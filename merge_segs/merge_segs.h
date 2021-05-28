#include "external.h"


void* ReadAllPixelsFromSingleBandRaster(string strRaster, unsigned int &nWidth, unsigned int &nHeight)
{

	GDALDataset* poInputDS = (GDALDataset*)GDALOpen((strRaster).c_str(), GA_ReadOnly);
	if (!poInputDS) return 0;
	//poInputDS->GetRasterXSize()
	//poInputDS->GetGCPSpatialRef();

	void* pabOutput = 0;

	nWidth = poInputDS->GetRasterXSize();
	nHeight = poInputDS->GetRasterYSize();

	switch (poInputDS->GetRasterBand(1)->GetRasterDataType())
	{
	case GDT_Float32:
		pabOutput = new float[nWidth*nHeight];
		poInputDS->RasterIO(GF_Read, 0, 0, nWidth, nHeight, pabOutput, nWidth, nHeight, GDT_Float32, 1, 0, 0, 0, 0, 0);
		break;
	case GDT_UInt32:
		pabOutput = new unsigned int[nWidth*nHeight];
		poInputDS->RasterIO(GF_Read, 0, 0, nWidth, nHeight, pabOutput, nWidth, nHeight, GDT_UInt32, 1, 0, 0, 0, 0, 0);
		break;
	case GDT_Int32:
		pabOutput = new unsigned int[nWidth*nHeight];
		poInputDS->RasterIO(GF_Read, 0, 0, nWidth, nHeight, pabOutput, nWidth, nHeight, GDT_UInt32, 1, 0, 0, 0, 0, 0);
		break;
	default:
		pabOutput = new uint16_t[nWidth*nHeight];
		poInputDS->RasterIO(GF_Read, 0, 0, nWidth, nHeight, pabOutput, nWidth, nHeight, GDT_UInt16, 1, 0, 0, 0, 0, 0);
		break;
	}


	GDALClose(poInputDS);
	return pabOutput;
}

GDALDataset* CreateGDALDatasetForWritingUsingPattern(string strNewRaster, string strExistingRaster)
{

	GDALDataset* poInputDS = (GDALDataset*)GDALOpen((strExistingRaster).c_str(), GA_ReadOnly);
	if (!poInputDS) return 0;


	void* pabOutput = 0;

	unsigned int nWidth = poInputDS->GetRasterXSize();
	unsigned int nHeight = poInputDS->GetRasterYSize();
	double padfGeoTransform[6];
	poInputDS->GetGeoTransform(padfGeoTransform);

	char **papszOptions = NULL;
	papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");



	GDALDataset* poOutputDS = (GDALDataset*)GDALCreate(
		GDALGetDriverByName("GTiff"),
		strNewRaster.c_str(),
		poInputDS->GetRasterXSize(),
		poInputDS->GetRasterYSize(),
		1, poInputDS->GetRasterBand(1)->GetRasterDataType(), papszOptions);



	const char* strProjRef = poInputDS->GetProjectionRef();

	poOutputDS->SetProjection(strProjRef);

	//poOutputDS->SetSpatialRef(poInputDS->GetSpatialRef());
	poOutputDS->SetGeoTransform(padfGeoTransform);
	poOutputDS->GetRasterBand(1)->SetNoDataValue(0);

	GDALClose(poInputDS);
	CSLDestroy(papszOptions);

	return poOutputDS;
}



struct Pixel
{
	unsigned int m_nX;
	unsigned int m_nY;
	float m_dblH;
	Pixel(unsigned int nX, unsigned int nY, float dblH)
	{
		m_nX = nX;
		m_nY = nY;
		m_dblH = dblH;
	};
};

struct BorderSegment
{
	//unsigned int nNeighbourID;
	set<Pixel*> setPixels;
	bool Contains(Pixel* poPix)
	{
		return setPixels.find(poPix) != setPixels.end();
	};

	bool Add(Pixel* poPix)
	{
		if (Contains(poPix)) return false;
		setPixels.insert(poPix);
		return true;
	};

	unsigned int Size()
	{
		return setPixels.size();
	};
	//set
};

struct SegmentMeta
{
	unsigned int nID;
	unsigned int nArea;
	unsigned int nCropPixelsCount;
	unsigned int nSampleX;
	unsigned int nSampleY;
	unsigned int nBorderLength;
	std::map<unsigned int, BorderSegment> mapBorderSegments;
	unsigned int nMinX, nMinY, nMaxX, nMaxY;


	bool IsRoughlyCrop(bool bStrong = false)
	{
		if (this->nID == 1) return false; //we should drop the marginal segment



		if (nArea == 0) return false;

		float dblElongation = ((float)nBorderLength) / nArea;


		if (!bStrong)
		{
			if (nArea == 0) return false;
			if (nArea < 10) return ((float)nCropPixelsCount) / nArea > 0.7999 ? true : false;
			if (nArea < 50) return ((float)nCropPixelsCount) / nArea > 0.71999 ? true : false;
			else if (dblElongation > 0.35) return ((float)nCropPixelsCount) / nArea > 0.71999 ? true : false;
			else
			{
				if (nArea < 250) return ((float)nCropPixelsCount) / nArea > 0.66 ? true : false;
				else return ((float)nCropPixelsCount) / nArea > 0.55 ? true : false;
			}
		}
		else
		{
			return ((float)nCropPixelsCount) / nArea > 0.9 ? true : false;
		}

	};

	SegmentMeta()
	{
		nArea = 0;
		nID = 0;
		nCropPixelsCount = 0;
		nSampleX = (nSampleY = 0);
		nBorderLength = 0;
		nMinX = (nMinY = 1e+6);
		nMaxX = (nMaxY = 0);

	};

	struct CompareSegments {
		bool operator()(const SegmentMeta* a, const SegmentMeta* b) const {

			return a->nArea < b->nArea ? true :
				a->nArea == b->nArea ? a->nID < b->nID : false;
		}
	};

};

class SegmentCatalogue
{
public:

	bool DropAllNoCrop()
	{
		list<SegmentMeta*> listSegToDrop;

		for (auto el : m_mapSegments)
		{
			if (!el.second->IsRoughlyCrop()) listSegToDrop.push_back(el.second);
		}

		SegmentMeta* poNoCropSeg = new SegmentMeta();
		AddSegment(poNoCropSeg);


		for (auto el : listSegToDrop)
		{
			MergePair(el, poNoCropSeg);
		}
		return true;
	};

	/*
	vector<unsigned int> FindLargeNeighbours(SegmentMeta *poSeg,
												unsigned int nMina)
	{
		unsigned int nSize = 25;
		SegmentMeta* poNeighbour;
		vector<unsigned int> vecLargeNeighbours;

		for (auto el : poSeg->mapNeighbours)
		{
			poNeighbour = GetSegmentRef(el.first);
			if (poNeighbour->nArea >= nMina && poNeighbour->IsRoughlyCrop())
				vecLargeNeighbours.push_back(el.first);
		}

		return vecLargeNeighbours;
	};
	*/
	/*
	bool ProcessBetweenTwoLargeCase(SegmentMeta* poSeg, unsigned int nMina)
	{
		return false;
		vector<unsigned int> vecLargeNeighbours = FindLargeNeighbours(poSeg, nMina);

		if (vecLargeNeighbours.size() != 2)
		{
			vecLargeNeighbours.clear();
			return false;
		}



		float dblInñl = (poSeg->mapNeighbours[vecLargeNeighbours[0]] + poSeg->mapNeighbours[vecLargeNeighbours[1]])
			/ ((float)poSeg->nBorderLength);

		//if (poSeg->nArea >= 25) return dblInñl >= 0.75 ? true : false;
		//else
		if (dblInñl < 0.99) return false;
		else
		{
			if (poSeg->mapNeighbours[vecLargeNeighbours[0]] > poSeg->mapNeighbours[vecLargeNeighbours[1]])
				MergePair(poSeg, GetSegmentRef(vecLargeNeighbours[0]));
			else
				MergePair(poSeg, GetSegmentRef(vecLargeNeighbours[1]));
			return true;
		}
	};
	*/




	//new merge strategy:
	// loop all segments which area greater than mina
	// try merge each neighbour if:
	// - border height between them less than threshold
	// - common border / all border less than
	// - segment inside convex hull or inside
	unsigned int FilterNotMergedSmallSegments(unsigned int nMinArea)
	{
		//we want to sellect small crop segments that satisfies (or):
		// - borders with two large segment and less than 25 - merge to one of them
		// - borders (>0.75) with two large segment and larger than or equal to 25 -  don't merge, don't delete
		// - group of segments wich can be merged together and then kept if and only if:
		//     inclusiveness in two large segments 0.75)
		// - 

		list<SegmentMeta*> listSegsToProcess;
		for (auto el : m_mapSegments)
		{
			if (el.second->nArea < nMinArea)
				listSegsToProcess.push_back(el.second);
		}

		list<SegmentMeta*> listSegsToDelete;
		for (auto el : listSegsToProcess)
			//if (!ProcessBetweenTwoLargeCase(el, nMinArea))
			listSegsToDelete.push_back(el);

		for (auto el : listSegsToDelete)
			RemoveSegment(el);
		/*
		SegmentMeta* poSegMeta = GetSmallest();
		unsigned int nCount = 0;
		while (poSegMeta->nArea < nArea && poSegMeta)
		{
			if (RemoveSegment(poSegMeta))
				nCount++;
			poSegMeta = GetSmallest();
		}
		*/

		//return nCount;
		return 0;
	}

	bool MergeIfInsideCropGroup(unsigned int nMinArea)
	{

		std::set<unsigned int> setSmallSegs;
		for (auto el : this->m_mapSegments)
			if (el.second->nArea < nMinArea) setSmallSegs.insert(el.first);

		for (auto nID : setSmallSegs)
		{
			SegmentMeta* poSeg = GetSegmentRef(nID);

			SegmentMeta* poLargeNeighbour = 0;
			bool bWasNoneCrop = false;
			bool bDontMerge = false;
			for (auto el : poSeg->mapBorderSegments)
			{
				SegmentMeta* poNeighbour = GetSegmentRef(el.first);
				if (!poNeighbour) continue;
				else if (!poNeighbour->IsRoughlyCrop())
				{
					if ((poSeg->nBorderLength < 10) || bWasNoneCrop || (el.second.Size() > 2))
					{
						bDontMerge = true;
						break;
					}
					else
						bWasNoneCrop = true;
				}
				else if (poNeighbour->nArea >= nMinArea)
				{
					if (poLargeNeighbour)
					{
						bDontMerge = true;
						break;
					}
					else poLargeNeighbour = poNeighbour;
				}
				//merge segment small segment to a large one if:
				//all neighbours are crop except one or two pixels
				//only one large crop neighbour
			}
			if ((!bDontMerge) && poLargeNeighbour)
			{
				MergePair(poSeg, poLargeNeighbour);
				//std::cout << poSeg->nID << std::endl;
			}
		}

		return true;

	};


	bool TryToBuildUpSegment(unsigned int nID,
		int nMinArea,
		float dblWeakHeight,
		float dblStrongHeight)
	{
		SegmentMeta* poSeg = GetSegmentRef(nID);
		SegmentMeta* poNeighbourSeg;
		float dblH, dblWeakIncl, dblStrongIncl;
		list<SegmentMeta*> listNeighrboursToMerge;

		while (true)
		{
			for (auto el : poSeg->mapBorderSegments)
			{
				poNeighbourSeg = GetSegmentRef(el.first);

				//if (poNeighbourSeg->nID == 0 || poNeighbourSeg->nArea >= nMinArea) continue;
				if (poNeighbourSeg->nID == 0 || poNeighbourSeg->nArea > poSeg->nArea) continue;

				bool bMerge = false;
				if (IsOneInsideAnother(poSeg, poNeighbourSeg))
					bMerge = true;
				//else
				//{
				dblH = CalculateBorderHeightBetweenSegments(poSeg, poNeighbourSeg);
				dblWeakIncl = CalcWeakInclusiveness(poSeg, poNeighbourSeg);
				dblStrongIncl = CalcStrongInclusiveness(poSeg, poNeighbourSeg);


				if (dblStrongIncl > 0.5001)
				{
					bMerge = true;
				}
				else if (dblH <= dblStrongHeight)
				{
					if (dblH < dblWeakHeight)
					{
						if ((dblWeakIncl >= 0.66)
							|| ((dblStrongIncl > 0.33) && (poNeighbourSeg->IsRoughlyCrop(true)))
							)
							bMerge = true;
					}
					else
					{
						if ((dblWeakIncl > 0.5)
							&& (dblStrongIncl > 0.33) && (poNeighbourSeg->IsRoughlyCrop(true)))
							bMerge = true;
					}
				}
				//}

				if (bMerge) listNeighrboursToMerge.push_back(poNeighbourSeg);

			}

			if (listNeighrboursToMerge.size() == 0) break;
			else
			{
				for (auto el : listNeighrboursToMerge)
				{
					MergePair(el, poSeg);
				}
				listNeighrboursToMerge.clear();
			}
		}

		return true;
	};

	bool RunBatchMerging(int nMinArea,
		float dblWeakHeight,
		float dblStrongHeight)
	{

		list<unsigned int> listSegsToProcess;
		for (auto el : this->m_mapSegments)
		{
			if (el.second->nArea >= nMinArea && el.second->IsRoughlyCrop())
				listSegsToProcess.push_back(el.second->nID);
		}

		for (auto el : listSegsToProcess)
		{
			if (this->GetSegmentRef(el))
				TryToBuildUpSegment(el, nMinArea, dblWeakHeight, dblStrongHeight);
		}

		return true;
	};

	unsigned int GetAfterMergeSegmentID(unsigned int nID)
	{

		std::map<unsigned int, unsigned int>::iterator iter;
		while ((iter = m_mapMergedSegments.find(nID)) != m_mapMergedSegments.end())
			nID = (*iter).second;


		return this->GetSegmentRef(nID) != 0 ? nID : 0;
	}

	void EraseFromMergedSegments(unsigned int nID)
	{
		m_mapMergedSegments.erase(nID);
	}


	bool SaveMergedSegments(string strOutputRaster, string strInputRaster, bool bDeleteMarginal = false, unsigned int nAddVal = 0)
	{
		GDALDataset* poInputDS = (GDALDataset*)GDALOpen((strInputRaster).c_str(), GA_ReadOnly);
		//poInputDS->GetRasterXSize()
		//poInputDS->GetGCPSpatialRef();
		double dblGeotransform[6];
		poInputDS->GetGeoTransform(dblGeotransform);

		if (MPLFileSys::FileExists(strOutputRaster))
			MPLFileSys::FileDelete(strOutputRaster);

		char **papszOptions = NULL;
		papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");

		GDALDataset* poOutputDS = (GDALDataset*)GDALCreate(
			GDALGetDriverByName("GTiff"),
			strOutputRaster.c_str(),
			poInputDS->GetRasterXSize(),
			poInputDS->GetRasterYSize(),
			1, GDT_UInt32, papszOptions);

		/*
		const char* strProjRef      = this->p_gdal_ds_->GetProjectionRef();

	if (OGRERR_NONE == srs.SetFromUserInput(strProjRef)) return true;
	else if (MPLFileSys::FileExists(MPLFileSys::RemoveExtension(this->raster_file_)+".prj"))
	{
		string prjFile		= MPLFileSys::RemoveExtension(this->raster_file_)+".prj";
		if (OGRERR_NONE==srs.SetFromUserInput(prjFile.c_str())) return true;
	}
	else if (MPLFileSys::FileExists(MPLFileSys::RemoveExtension(this->raster_file_)+".tab"))
	{
		string tabFile = MPLFileSys::RemoveExtension(this->raster_file_)+".tab";
		if (ReadSpatialRefFromMapinfoTabFile(tabFile,srs)) return true;
	}
		*/

		const char* strProjRef = poInputDS->GetProjectionRef();
		//OGRSpatialReference oSRS;
		//oSRS.SetFromUserInput(strProjRef);
		poOutputDS->SetProjection(strProjRef);

		//poOutputDS->SetSpatialRef(poInputDS->GetSpatialRef());
		poOutputDS->SetGeoTransform(dblGeotransform);
		poOutputDS->GetRasterBand(1)->SetNoDataValue(0);

		//poOutputDS->set

		int nChunkSize = 256;
		for (int nPosY = 0; nPosY < poInputDS->GetRasterYSize(); nPosY += nChunkSize)
		{
			int nChunkHeight = nPosY + nChunkSize < poInputDS->GetRasterYSize() ? nChunkSize :
				poInputDS->GetRasterYSize() - nPosY;

			int nArea = nChunkHeight * poInputDS->GetRasterXSize();
			unsigned int* panValues = new unsigned int[nArea];
			poInputDS->RasterIO(GF_Read, 0, nPosY, poInputDS->GetRasterXSize(), nChunkHeight,
				panValues, poInputDS->GetRasterXSize(), nChunkHeight,
				GDT_UInt32, 1, 0, 0, 0, 0);

			unsigned int nNewVal;
			SegmentMeta* poSeg;
			for (int i = 0; i < nArea; i++)
			{
				panValues[i] = GetSegmentRef(panValues[i]) ? panValues[i] : GetAfterMergeSegmentID(panValues[i]);

				if (bDeleteMarginal)
				{
					if (panValues[i])
					{
						poSeg = GetSegmentRef(panValues[i]);
						if ((poSeg->nMinX == 0) || (poSeg->nMaxX == poInputDS->GetRasterXSize() - 1)
							|| (poSeg->nMinY == 0) || (poSeg->nMaxY == poInputDS->GetRasterYSize() - 1))
						{
							panValues[i] = 0;
						}
					}
				}

				if (panValues[i]) panValues[i] += nAddVal;

			}

			poOutputDS->RasterIO(GF_Write, 0, nPosY, poInputDS->GetRasterXSize(), nChunkHeight,
				panValues, poInputDS->GetRasterXSize(), nChunkHeight,
				GDT_UInt32, 1, 0, 0, 0, 0);
			delete(panValues);

		}

		GDALClose(poInputDS);
		GDALClose(poOutputDS);
		CSLDestroy(papszOptions);

		return true;
	}

	bool InitFromFiles(string strSegmentedImage, string strCropMaskImage, string strHeightsImage = "")
	{
		const int nChunkMaxHeight = 256;

		GDALDataset* poSegmentedDS = (GDALDataset*)GDALOpen((strSegmentedImage).c_str(), GA_ReadOnly);
		if (!poSegmentedDS) return false;

		GDALDataset* poCropMaskDS = (GDALDataset*)GDALOpen((strCropMaskImage).c_str(), GA_ReadOnly);
		if (!poCropMaskDS) return false;

		GDALDataset* poHeightsDS = (strHeightsImage == "") ? 0
			: (GDALDataset*)GDALOpen((strHeightsImage).c_str(), GA_ReadOnly);

		for (int nPosY = 0; nPosY < poSegmentedDS->GetRasterYSize(); nPosY += nChunkMaxHeight)
		{
			ProcessDataChunk(nPosY, nChunkMaxHeight, poSegmentedDS, poCropMaskDS, poHeightsDS);
		}

		GDALClose(poSegmentedDS);
		GDALClose(poCropMaskDS);
		return true;
	};

	bool RemoveSegment(SegmentMeta* poSegMeta)
	{
		SegmentMeta* poNeighborSeg;
		for (auto el : poSegMeta->mapBorderSegments)
		{
			if (poNeighborSeg = GetSegmentRef(el.first))
			{
				poNeighborSeg->mapBorderSegments.erase(el.first);
			}
		}
		/*
		for (auto el : poSegMeta->mapNeighbours)
		{
			if (poNeighborSeg = GetSegmentRef(el.first))
			{
				if (poNeighborSeg->mapNeighbours.find(poSegMeta->nID) != poNeighborSeg->mapNeighbours.end())
					poNeighborSeg->mapNeighbours.erase(poSegMeta->nID);
				if (poNeighborSeg->mapBorderHeightsByNeighbour.find(poSegMeta->nID)
						!= poNeighborSeg->mapBorderHeightsByNeighbour.end())
					poNeighborSeg->mapBorderHeightsByNeighbour.erase(poSegMeta->nID);
			}
		}
		*/


		m_mapSegments.erase(poSegMeta->nID);
		m_setAreaIndex.erase(poSegMeta);
		delete(poSegMeta);

		return true;
	};

	SegmentMeta* GetSegmentRef(unsigned int nID)
	{
		std::map<unsigned int, SegmentMeta*>::iterator iter;
		return (iter = m_mapSegments.find(nID)) == m_mapSegments.end() ? 0 : (*iter).second;
	};

	bool AddSegment(SegmentMeta* poSegMeta)
	{
		m_mapSegments[poSegMeta->nID] = poSegMeta;
		m_setAreaIndex.insert(poSegMeta);
		return true;
	};

	unsigned int GetSmallestID()
	{
		return m_mapSegments.size() > 0 ? (*m_setAreaIndex.begin())->nID : 0;
	};

	SegmentMeta* GetSmallest()
	{
		return m_mapSegments.size() > 0 ? (*m_setAreaIndex.begin()) : 0;
	};

	unsigned int GetLargestID()
	{
		return m_mapSegments.size() > 0 ? (*m_setAreaIndex.rbegin())->nID : 0;
	};

	unsigned int GetSize()
	{
		return GetSegmentRef(0) != 0 ? m_mapSegments.size() - 1 : m_mapSegments.size();
	};

	bool UpdateSegmentArea(const unsigned int &nID, const unsigned int &nArea)
	{
		SegmentMeta* poSeg = m_mapSegments[nID];
		m_setAreaIndex.erase(poSeg);
		poSeg->nArea = nArea;
		m_setAreaIndex.insert(poSeg);
		return true;
	};

	~SegmentCatalogue()
	{
		for (auto seg : m_mapSegments)
			delete(seg.second);
	}

protected:
	bool IsBboxInside(SegmentMeta* poExternalSeg, OGREnvelope &oEnvp)
	{
		double E = 1e-4;
		return ((oEnvp.MinX > poExternalSeg->nMinX + E) &&
			(oEnvp.MinY > poExternalSeg->nMinY + E) &&
			(oEnvp.MaxX < poExternalSeg->nMaxX - E) &&
			(oEnvp.MaxY < poExternalSeg->nMaxY - E)) ? true : false;
	}


	bool IsOneInsideAnother(SegmentMeta* poExternalSeg, SegmentMeta* poInternalSeg)
	{
		std::set<unsigned int> setProcessedSegs;
		return IsOneInsideAnotherRecursive(poExternalSeg, poInternalSeg, &setProcessedSegs);
	}


	bool IsOneInsideAnotherRecursive(SegmentMeta* poExternalSeg,
		SegmentMeta* poInternalSeg,
		std::set<unsigned int>* posetProcessedSegs)
	{

		OGREnvelope oEnvp;
		oEnvp.MinX = poInternalSeg->nMinX;
		oEnvp.MinY = poInternalSeg->nMinY;
		oEnvp.MaxX = poInternalSeg->nMaxX;
		oEnvp.MaxY = poInternalSeg->nMaxY;
		if (!IsBboxInside(poExternalSeg, oEnvp))
			return false;


		for (auto el : poInternalSeg->mapBorderSegments)
		{
			if (el.first == poExternalSeg->nID) continue;
			else if (posetProcessedSegs->find(el.first) != posetProcessedSegs->end()) continue;
			else
			{
				posetProcessedSegs->insert(poInternalSeg->nID);
				if (!IsOneInsideAnotherRecursive(poExternalSeg, GetSegmentRef(el.first), posetProcessedSegs))
					return false;
			}
		}

		return true;
	}

	float CalcWeakInclusiveness(SegmentMeta* poSeg1, SegmentMeta* poSeg2)
	{
		return poSeg2->mapBorderSegments.find(poSeg1->nID) == poSeg2->mapBorderSegments.end() ? 0
			: ((float)poSeg2->mapBorderSegments[poSeg1->nID].Size()) / poSeg2->nBorderLength;
		/*
		return poSeg2->mapNeighbours.find(poSeg1->nID) != poSeg2->mapNeighbours.end() ?
			((float)poSeg2->mapNeighbours[poSeg1->nID])/poSeg2->nBorderLength: 0;
		*/
	}

	float CalcStrongInclusiveness(SegmentMeta* poSeg1, SegmentMeta* poSeg2)
	{
		if (poSeg2->mapBorderSegments.find(poSeg1->nID) == poSeg2->mapBorderSegments.end())
			return 0;

		unsigned int nExclusiveBoundaryPixelsCount = 0;
		bool bExclusive;

		for (auto poPix : poSeg2->mapBorderSegments[poSeg1->nID].setPixels)
		{
			bExclusive = true;
			for (auto el : poSeg2->mapBorderSegments)
			{
				if (el.first == poSeg1->nID) continue;
				if (el.second.Contains(poPix))
				{
					bExclusive = false;
					break;
				}
			}
			if (bExclusive) nExclusiveBoundaryPixelsCount++;
		}
		/*
		for (auto el : poSeg2->mapBorderSegments)
		{
			if (el.first == poSeg1->nID) continue;

			for (auto poPix : el.second.setPixels)
			{
				if (!poSeg2->mapBorderSegments[poSeg1->nID].Contains(poPix))
					nOtherNeighboursPixCount++;
			}
		}
		*/

		return ((float)nExclusiveBoundaryPixelsCount) / poSeg2->nBorderLength;
		//	((float)poSeg2->mapBorderSegments[poSeg1->nID].Size()) /
		//	(poSeg2->mapBorderSegments[poSeg1->nID].Size() + nOtherNeighboursPixCount);

		/*
		if (poSeg2->mapNeighbours.find(poSeg1->nID) == poSeg2->mapNeighbours.end()) return 0;

		unsigned int nOtherNeighboursLen = 0;
		for (auto el : poSeg2->mapNeighbours)
			if (el.first != poSeg1->nID) nOtherNeighboursLen += el.second;

		return ((float)poSeg2->mapNeighbours[poSeg1->nID])
				/(poSeg2->mapNeighbours[poSeg1->nID] + nOtherNeighboursLen);
				*/
	}

	float CalculateBorderHeightBetweenSegments(SegmentMeta* poSeg1, SegmentMeta* poSeg2)
	{
		if (poSeg1->mapBorderSegments.find(poSeg2->nID) == poSeg1->mapBorderSegments.end())
			return 0;

		list<float> listHeights;
		for (auto poPix : poSeg1->mapBorderSegments[poSeg2->nID].setPixels)
			listHeights.push_back(poPix->m_dblH);

		for (auto poPix : poSeg2->mapBorderSegments[poSeg1->nID].setPixels)
			listHeights.push_back(poPix->m_dblH);

		listHeights.sort();

		std::list<float>::iterator it = listHeights.begin();
		std::advance(it, (listHeights.size() + 1) / 2);

		return (*it);


	}

	/*
	bool TryToMergeSegment(SegmentMeta* poSeg, float dblMaxHeight = 0)
	{
		//if (!poSeg->IsRoughlyCrop()) return false;

		unsigned int nMaxLength = 0;
		float dblMinNeighbourHeight = 1.1;
		float dblH;
		unsigned int nMaxNeighbourLength = 0;
		unsigned int nBestNeighbour = 0;
		SegmentMeta* poSegNeighbour = 0;
		for (auto el : poSeg->mapBorderSegments)
		{
			poSegNeighbour = GetSegmentRef(el.first);
			if (poSegNeighbour->IsRoughlyCrop())
			{
				if (dblMaxHeight > 0)
				{
					dblH = CalculateBorderHeightBetweenSegments(poSeg, poSegNeighbour);
					if ((dblH > 0) && (dblH < dblMinNeighbourHeight))
					{
						dblMinNeighbourHeight = dblH;
						nBestNeighbour = el.first;
					}
				}
				else
				{
					if (el.second > nMaxNeighbourLength)
					{
						nMaxNeighbourLength = el.second;
						nBestNeighbour = el.first;
					}
				}
			}
		}

		//std::cout << dblMinNeighbourHeight << " " << dblMaxHeight << endl;
		if (dblMaxHeight > 0)
			return dblMinNeighbourHeight < dblMaxHeight ? MergePair(poSeg, GetSegmentRef(nBestNeighbour)) : false;
		else
			return nMaxNeighbourLength > 0 ? MergePair(poSeg, GetSegmentRef(nBestNeighbour)) : false;

	};
	*/

	bool MergePair(SegmentMeta* poSegRemoved, SegmentMeta* poSegMerged)
	{
		//TODO: merge  mapBorderHeightsByNeighbour
		//
		//std::cout << "MergePair" << endl;

		unsigned int nIDRemoved = poSegRemoved->nID;
		unsigned int nIDMerged = poSegMerged->nID;


		m_mapSegments.erase(poSegRemoved->nID);
		m_setAreaIndex.erase(poSegRemoved);
		m_setAreaIndex.erase(poSegMerged);

		poSegMerged->nArea += poSegRemoved->nArea;
		poSegMerged->nCropPixelsCount += poSegRemoved->nCropPixelsCount;
		if (poSegMerged->mapBorderSegments.find(nIDRemoved) != poSegMerged->mapBorderSegments.end())
			poSegMerged->mapBorderSegments.erase(nIDRemoved);
		//poSegMerged->mapNeighbours.erase(nIDRemoved);
		//poSegMerged->mapBorderHeightsByNeighbour.erase(nIDRemoved);

		SegmentMeta* poSegNeighbour = 0;
		for (auto el : poSegRemoved->mapBorderSegments)
		{
			if (el.first == nIDMerged) continue;

			poSegNeighbour = GetSegmentRef(el.first);

			if (poSegMerged->mapBorderSegments.find(el.first) != poSegMerged->mapBorderSegments.end())
			{
				for (auto poPix : poSegRemoved->mapBorderSegments[el.first].setPixels)
					poSegMerged->mapBorderSegments[el.first].Add(poPix);

				if (el.first != 0)
				{
					for (auto poPix : poSegNeighbour->mapBorderSegments[nIDRemoved].setPixels)
						poSegNeighbour->mapBorderSegments[nIDMerged].Add(poPix);
				}
			}
			else
			{
				poSegMerged->mapBorderSegments[el.first] = el.second;
				if (el.first != 0)
				{
					poSegNeighbour->mapBorderSegments[nIDMerged]
						= poSegNeighbour->mapBorderSegments[nIDRemoved];
				}
			}
			if (el.first > 0) poSegNeighbour->mapBorderSegments.erase(nIDRemoved);


		}

		if (poSegMerged->nMaxX < poSegRemoved->nMaxX) poSegMerged->nMaxX = poSegRemoved->nMaxX;
		if (poSegMerged->nMaxY < poSegRemoved->nMaxY) poSegMerged->nMaxY = poSegRemoved->nMaxY;
		if (poSegMerged->nMinX > poSegRemoved->nMinX) poSegMerged->nMinX = poSegRemoved->nMinX;
		if (poSegMerged->nMinY > poSegRemoved->nMinY) poSegMerged->nMinY = poSegRemoved->nMinY;


		m_setAreaIndex.insert(poSegMerged);

		delete(poSegRemoved);

		m_mapMergedSegments[nIDRemoved] = nIDMerged;
		return true;

	};
	bool ProcessDataChunk(int nPosY,
		int nChunkMaxHeight,
		GDALDataset* poSegmentedDS,
		GDALDataset* poCropMaskDS,
		GDALDataset* poHeightsDS = 0)
	{
		int nWidth = poSegmentedDS->GetRasterXSize();

		int nChunkHeight = nPosY + nChunkMaxHeight > poSegmentedDS->GetRasterYSize() ?
			poSegmentedDS->GetRasterYSize() - nPosY : nChunkMaxHeight;

		nChunkHeight += 2;
		int nChunkArea = nWidth * nChunkHeight;
		unsigned int* panPixels = new unsigned int[nChunkArea];
		uint16_t* panCropPixels = new uint16_t[nChunkArea];
		float* padblHeights = (poHeightsDS == 0) ? 0 : new float[nChunkArea];


		for (int i = 0; i < nChunkArea; i++)
			panPixels[i] = 0;


		bool bFirstChunk = nPosY == 0 ? true : false;
		bool bEdgeChunk = bFirstChunk ? true :
			(nPosY + nChunkHeight - 2) == poSegmentedDS->GetRasterYSize() ? true : false;

		poSegmentedDS->RasterIO(GF_Read, 0,
			bFirstChunk ? 0 : nPosY - 1, poSegmentedDS->GetRasterXSize(),
			bEdgeChunk ? nChunkHeight - 1 : nChunkHeight,
			bFirstChunk ? &panPixels[nWidth] : panPixels,
			poSegmentedDS->GetRasterXSize(),
			bEdgeChunk ? nChunkHeight - 1 : nChunkHeight,
			GDT_UInt32, 1, 0, 0, 0, 0);
		poCropMaskDS->RasterIO(GF_Read, 0,
			bFirstChunk ? 0 : nPosY - 1, poSegmentedDS->GetRasterXSize(),
			bEdgeChunk ? nChunkHeight - 1 : nChunkHeight,
			bFirstChunk ? &panCropPixels[nWidth] : panCropPixels,
			poSegmentedDS->GetRasterXSize(),
			bEdgeChunk ? nChunkHeight - 1 : nChunkHeight,
			GDT_UInt16, 1, 0, 0, 0, 0);
		if (poHeightsDS)
		{
			poHeightsDS->RasterIO(GF_Read, 0,
				bFirstChunk ? 0 : nPosY - 1, poSegmentedDS->GetRasterXSize(),
				bEdgeChunk ? nChunkHeight - 1 : nChunkHeight,
				bFirstChunk ? &padblHeights[nWidth] : padblHeights,
				poSegmentedDS->GetRasterXSize(),
				bEdgeChunk ? nChunkHeight - 1 : nChunkHeight,
				GDT_Float32, 1, 0, 0, 0, 0);
		}

		unsigned int panNeighbours[4];
		for (int i = 1; i < nChunkHeight - 1; i++)
		{
			for (int j = 0; j < nWidth; j++)
			{
				if (!panPixels[nWidth*i + j]) continue;

				panNeighbours[0] = panPixels[nWidth*(i - 1) + j];
				panNeighbours[2] = panPixels[nWidth*(i + 1) + j];
				panNeighbours[1] = j == nWidth - 1 ? 0 : panPixels[nWidth*i + j + 1];
				panNeighbours[3] = j == 0 ? 0 : panPixels[nWidth*i + j - 1];

				if (poHeightsDS == 0)
				{
					CollectPrimaryInfoInPixel(panPixels[nWidth*i + j], j, nPosY + i - 1,
						panCropPixels[nWidth*i + j], panNeighbours);
				}
				else
				{
					CollectAdditionalInfoInPixel(panPixels[nWidth*i + j], j, nPosY + i - 1,
						panNeighbours, padblHeights[nWidth*i + j]);
				}
			}
		}


		delete[]panPixels;
		delete[]panCropPixels;
		delete[]padblHeights;

		return true;
	};

	void CollectAdditionalInfoInPixel(const unsigned int &nID,
		const unsigned int nX,
		const unsigned int nY,
		const unsigned int panNeighbours[4],
		const float &dblHeight)
	{
		if (nID == 0 || GetAfterMergeSegmentID(nID) == 0) return;

		SegmentMeta* poSegMeta = GetSegmentRef(nID);

		Pixel *poPixel = 0;
		unsigned int nNeighbour;
		for (int i = 0; i < 4; i++)
		{
			if (panNeighbours[i] == 0) continue;

			if (panNeighbours[i] != poSegMeta->nID)
			{
				if (!poPixel)
					poPixel = new Pixel(nX, nY, dblHeight);

				nNeighbour = GetAfterMergeSegmentID(panNeighbours[i]) == 0 ? 0 : panNeighbours[i];


				if (poSegMeta->mapBorderSegments.find(nNeighbour) != poSegMeta->mapBorderSegments.end())
				{
					poSegMeta->mapBorderSegments[nNeighbour].Add(poPixel);
				}
				else
					(poSegMeta->mapBorderSegments[nNeighbour] = BorderSegment()).Add(poPixel);

			}
		}


	}

	void CollectPrimaryInfoInPixel(const unsigned int &nID,
		const int &nX,
		const int &nY,
		const int &nCropStatus,
		const unsigned int panNeighbours[4])
	{

		SegmentMeta* poSegMeta = GetSegmentRef(nID);
		if (!poSegMeta)
		{
			poSegMeta = new SegmentMeta();
			poSegMeta->nSampleX = nX;
			poSegMeta->nSampleY = nY;
			poSegMeta->nArea = 1;
			poSegMeta->nID = nID;
			AddSegment(poSegMeta);
		}
		else
			UpdateSegmentArea(nID, poSegMeta->nArea + 1);

		poSegMeta->nCropPixelsCount += nCropStatus;

		if (nX > poSegMeta->nMaxX) poSegMeta->nMaxX = nX;
		if (nY > poSegMeta->nMaxY) poSegMeta->nMaxY = nY;
		if (nX < poSegMeta->nMinX) poSegMeta->nMinX = nX;
		if (nY < poSegMeta->nMinY) poSegMeta->nMinY = nY;


		for (int i = 0; i < 4; i++)
		{
			if ((panNeighbours[i] != 0) && (panNeighbours[i] != poSegMeta->nID))
			{
				poSegMeta->nBorderLength++;
				return;
			}
		}

	};

protected:
	std::map<unsigned int, SegmentMeta*> m_mapSegments;
	std::set<SegmentMeta*, SegmentMeta::CompareSegments> m_setAreaIndex;
	std::map<unsigned int, unsigned int> m_mapMergedSegments;


};