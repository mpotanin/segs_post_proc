// mix_segmentations.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../merge_segs/merge_segs.h"

//init SegmentCatalogue with subsegment segmentation
//init SegmentationMixerCatalogue
//

struct SubSegmentMeta
{
	unsigned int m_nArea;
	unsigned int m_nID;
	//map<unsigned int, unsigned int> m_mapNeighbours;
	map<unsigned int, unsigned int> m_mapSuperSegments;
	map<unsigned int, unsigned int> m_mapSubstitutes;
	unsigned int m_nProperSuperSegment;

	SubSegmentMeta()
	{
		m_nArea = 0;
		m_nID = 0;
		m_nProperSuperSegment = 0;
	};


	unsigned int GetProperSuperSegmentID()
	{
		if (m_nProperSuperSegment > 0) return m_nProperSuperSegment;

		unsigned int nMaxIntersection = 0;
		for (auto el : m_mapSuperSegments)
		{
			if (el.second > nMaxIntersection)
			{
				nMaxIntersection = el.second;
				m_nProperSuperSegment = el.first;
			}
		}
		return m_nProperSuperSegment;
	}

};

//debug
unsigned int nNoSubstituteCount = 0;
unsigned int nSubstituteCount = 0;
//end-debug

struct TempSegment
{
	unsigned int m_nMinY;
	unsigned int m_nSubSegID;
	unsigned int m_nSuperSegID;
	map<unsigned int, unsigned int> m_mapNeighbours;
	set<unsigned int> m_setPixels;

	bool UpdateNeighbours(unsigned int panNeighboursIndex[4],
		unsigned int* panSubSegPixels,
		unsigned int* panSuperSegPixels,
		map<unsigned int, SubSegmentMeta*> &mapSubSegments)
	{

		for (int n = 0; n < 4; n++)
		{
			if (panNeighboursIndex[n] == -1) continue;
			if (panSubSegPixels[panNeighboursIndex[n]] == m_nSubSegID) continue;
			if (panSuperSegPixels[panNeighboursIndex[n]] != m_nSuperSegID) continue;
			if (mapSubSegments[panSubSegPixels[panNeighboursIndex[n]]]->GetProperSuperSegmentID() !=
				panSuperSegPixels[panNeighboursIndex[n]]) continue;

			if (m_mapNeighbours.find(panSubSegPixels[panNeighboursIndex[n]]) == m_mapNeighbours.end())
				m_mapNeighbours[panSubSegPixels[panNeighboursIndex[n]]] = 1;
			else
				m_mapNeighbours[panSubSegPixels[panNeighboursIndex[n]]] += 1;
		}
		return true;
	}

	bool MergeSegment(TempSegment* poMergingPart)
	{
		m_nMinY = m_nMinY > poMergingPart->m_nMinY ? poMergingPart->m_nMinY : m_nMinY;

		for (auto pix : poMergingPart->m_setPixels)
			m_setPixels.insert(pix);

		for (auto el : poMergingPart->m_mapNeighbours)
		{
			if (m_mapNeighbours.find(el.first) == m_mapNeighbours.end())
				m_mapNeighbours[el.first] = el.second;
			else
				m_mapNeighbours[el.first] += el.second;
		}

		return true;
	};

	unsigned int CloseTempSegment(unsigned int* panOutputPixels, unsigned nMaxSegID)
	{
		unsigned int nIDSubstituteID;

		if (m_mapNeighbours.size() == 0)
		{
			//debug
			nNoSubstituteCount++;
			//end-debug

			nIDSubstituteID = nMaxSegID;
			nMaxSegID++;
		}
		else
		{
			//debug
			nSubstituteCount++;
			//end-debug

			unsigned int nLongestBorder = 0;

			for (auto el : m_mapNeighbours)
			{
				if (el.second > nLongestBorder)
				{
					nLongestBorder = el.second;
					nIDSubstituteID = el.first;
				}
			}
		}

		for (auto ind : m_setPixels)
		{
			panOutputPixels[ind] = nIDSubstituteID;
		}

		return nMaxSegID;
	};

};


class SegmentMixerCatalogue
{

public:
	bool InitFromArrays(unsigned int* panSubSegPixels,
		unsigned int* panSuperSegPixels,
		unsigned int nWidth,
		unsigned int nHeight)
	{
		unsigned int nTotalPixels = nWidth * nHeight;
		for (unsigned int i = 0; i < nTotalPixels; i++)
		{
			if (m_mapSubSegments.find(panSubSegPixels[i]) == m_mapSubSegments.end())
			{
				SubSegmentMeta* poSubSeg = new SubSegmentMeta();
				poSubSeg->m_nID = panSubSegPixels[i];
				m_mapSubSegments[panSubSegPixels[i]] = poSubSeg;
			}

			/*
			if (m_mapSuperSegments.find(panSuperSegPixels[i]) == m_mapSuperSegments.end())
			{
				m_mapSuperSegments[panSuperSegPixels[i]] = set<unsigned int>();
			}
			*/

			SubSegmentMeta* poSubSeg = m_mapSubSegments[panSubSegPixels[i]];
			poSubSeg->m_nArea += 1;
			if (poSubSeg->m_mapSuperSegments.find(panSuperSegPixels[i]) == poSubSeg->m_mapSuperSegments.end())
			{
				poSubSeg->m_mapSuperSegments[panSuperSegPixels[i]] = 1;
			}
			else poSubSeg->m_mapSuperSegments[panSuperSegPixels[i]] += 1;

			/*
			unsigned int panNeighboursID[4];
			panNeighboursID[0] = i < nWidth ? panSubSegPixels[i] : panSubSegPixels[i - nWidth];
			panNeighboursID[1] = i >= nTotalPixels - nWidth ? panSubSegPixels[i] : panSubSegPixels[i + nWidth];
			panNeighboursID[2] = i % nWidth == 0 ? panSubSegPixels[i] : panSubSegPixels[i - 1];
			panNeighboursID[3] = (i+1) % nWidth == 0 ? panSubSegPixels[i] : panSubSegPixels[i + 1];


			for (int n = 0; n < 4; n++)
			{
				if (panNeighboursID[n] != panSubSegPixels[i])
				{
					if (poSubSeg->m_mapNeighbours.find(panNeighboursID[n]) == poSubSeg->m_mapNeighbours.end())
						poSubSeg->m_mapNeighbours[panNeighboursID[n]] = 1;
					else
						poSubSeg->m_mapNeighbours[panNeighboursID[n]] += 1;
				}
			}
			*/

			/*
			if (m_mapSuperSegments[panSuperSegPixels[i]].find(panSubSegPixels[i]) ==
				m_mapSuperSegments[panSuperSegPixels[i]].end())
			{
				m_mapSuperSegments[panSuperSegPixels[i]].insert(panSubSegPixels[i]);
			}
			*/

			/*
			set<unsigned int>* posetSubSegIDs = m_mapSuperSegments[panSuperSegPixels[i]];
			if (posetSubSegIDs->find(panSubSegPixels[i]) == posetSubSegIDs->end())
				posetSubSegIDs->insert(panSubSegPixels[i]);
			*/
		}



		return true;
	};

	/*
	bool CalcSubstitutes()
	{
		std::cout << m_mapSubSegments.size() << endl;
		std::cout << m_mapSuperSegments.size() << endl;

		set<unsigned int> setProperSubSegs;
		set<unsigned int> setToBeFixedSubSegs;
		SubSegmentMeta* poSubSeg;
		unsigned int nCount = 0;
		unsigned int nCountNoSubSegments = 0;
		for (auto el : m_mapSuperSegments)
		{
			//debug
			nCount++;
			if (nCount % 10000 == 0) std::cout << nCount << endl;
			//end-debug

			setProperSubSegs.clear();
			setToBeFixedSubSegs.clear();

			for (auto nSubSegID : el.second)
			{
				poSubSeg = m_mapSubSegments[nSubSegID];
				if (poSubSeg->GetProperSuperSegmentID() == el.first)
				{
					setProperSubSegs.insert(nSubSegID);
				}
				else
				{
					setToBeFixedSubSegs.insert(nSubSegID);
				}
			}
			if ((setProperSubSegs.size() == 0) || (setToBeFixedSubSegs.size() == 0)) continue;

			for (auto nSubSegID : setToBeFixedSubSegs)
			{
				poSubSeg = m_mapSubSegments[nSubSegID];

				for (auto nSubSegIDProper : setProperSubSegs)
				{
					if (poSubSeg->m_mapNeighbours.find(nSubSegIDProper) != poSubSeg->m_mapNeighbours.end())
					{
						poSubSeg->m_mapSubstitutes[el.first] = nSubSegIDProper;
					}
				}
				//poSubSeg->m_mapSubstitutes[el.first] = (*setProperSubSegs.begin());
			}
		}

		//debug
		std::cout << nCountNoSubSegments << std::endl;
		//end-debug
		return true;
	};
	*/

	//when reading data save bordersegments for neighbours
	//don't cal calcproper func
	//when saving ouput
	//analyze pixel value for subseg if not proper
	//put into open new subsegment store all neighbours
	//when new subsegment is closed set 
	//get only proper from subsegs proper
	//choose among them neighbours

	/*
	bool CalcSubstitutes()
	{
		//loop through supersegments
		//select that must and mustn't be substitute
		//debug
		unsigned int nTotalNoSubSegs = 0;
		unsigned int nCount = 0;
		//end-dbug

		std::cout << m_mapSubSegments.size() << endl;
		std::cout << m_mapSuperSegments.size() << endl;

		SubSegmentMeta* poSubSegment = 0;
		set<unsigned int> setSubSegsFilt;

		for (auto el : m_mapSubSegments)
		{
			//debug
			nCount++;
			if (nCount % 10000 == 0) std::cout << nCount << endl;
			//end-debug

			poSubSegment = el.second;
			if (poSubSegment->m_mapSuperSegments.size() == 1) continue;

			unsigned int nProperSuperSeg = poSubSegment->GetProperSuperSegmentID();

			for (auto el2 : poSubSegment->m_mapSuperSegments)
			{
				if (el2.first != nProperSuperSeg)
				{
					nTotalNoSubSegs++;


					//set<unsigned int> setSubSegs = *m_mapSuperSegments[el2.first];

					setSubSegsFilt.clear();


					for (auto el3 : m_mapSuperSegments[el2.first])
					{
						nTotalNoSubSegs+=el3;
						continue;

						if (m_mapSubSegments[el3]->GetProperSuperSegmentID()==el2.first)
						{
							setSubSegsFilt.insert(el3);
						}
					}

					/*
					if (setSubSegsFilt.size() == 0)
					{
						//debug
						nTotalNoSubSegs++;
						//std::cout << "No subsegments: " << el2.first << std::endl;
						//end-debug
					}

					else
					{
						poSubSegment->m_mapSubstitutes[el.first] = 0;

						unsigned int nSubstituteID = 0;
						unsigned int nLongestBorder = 0;
						for (auto el3 : setSubSegsFilt)
						{
							if (poSubSegment->m_mapNeighbours.find(el3) != poSubSegment->m_mapNeighbours.end())
							{
								if (poSubSegment->m_mapNeighbours[el3] > nLongestBorder)
								{
									nLongestBorder = poSubSegment->m_mapNeighbours[el3];
									nSubstituteID = el3;
								}
							}
						}

						if (nSubstituteID == 0)
						{
							std::cout << "NO substitute: " << poSubSegment->m_nID << " " << el2.first << std::endl;
						}
						else
						{
							poSubSegment->m_mapSubstitutes[el2.first] = nSubstituteID;
						}
					}

				}
			}

		}

		//debug
		std::cout << nTotalNoSubSegs << endl;
		//end-debug
		return true;
	};
	*/
	bool SaveAdjustedSubsegments(unsigned int* panSubSegPixels,
		unsigned int* panSuperSegPixels,
		unsigned int nWidth,
		unsigned int nHeight,
		GDALDataset* poOutputDS)
	{
		unsigned int nTotalPixels = nWidth * nHeight;
		unsigned int* panOutputPixels = new unsigned int[nTotalPixels];
		SubSegmentMeta* poSubSeg = 0;
		set<TempSegment*> setTempSegs;
		unsigned int panNeighboursIndex[4];
		unsigned int nIncrementedMaxID = m_mapSubSegments.size() + 1;

		//debug
		std::cout << m_mapSubSegments.size() << std::endl;
		//end-debug

		for (int i = 0; i < nTotalPixels; i++)
		{
			poSubSeg = m_mapSubSegments[panSubSegPixels[i]];
			if (poSubSeg->GetProperSuperSegmentID() == panSuperSegPixels[i])
			{
				panOutputPixels[i] = panSubSegPixels[i];
			}
			else
			{
				panNeighboursIndex[0] = i < nWidth ? -1 : i - nWidth;
				panNeighboursIndex[1] = i % nWidth == 0 ? -1 : i - 1;
				panNeighboursIndex[2] = i + nWidth >= nTotalPixels ? -1 : i + nWidth;
				panNeighboursIndex[3] = (i + 1) % nWidth == 0 ? -1 : i + 1;

				TempSegment* paoTempSegs[2] = { 0,0 };
				for (int n = 0; n < 2; n++)
				{
					if (panNeighboursIndex[n] == -1) continue;
					if ((panSubSegPixels[panNeighboursIndex[n]] == panSubSegPixels[i]) &&
						(panSuperSegPixels[panNeighboursIndex[n]] == panSuperSegPixels[i]))
					{
						for (auto poTempSegIter : setTempSegs)
						{
							if ((poTempSegIter->m_nSubSegID == panSubSegPixels[panNeighboursIndex[n]]) &&
								(poTempSegIter->m_nSuperSegID == panSuperSegPixels[panNeighboursIndex[n]]))
							{
								if (poTempSegIter->m_setPixels.find(panNeighboursIndex[n]) !=
									poTempSegIter->m_setPixels.end())
									paoTempSegs[n] = poTempSegIter;
							}
						}
					}
				}

				TempSegment* poTempSeg = 0;
				if ((!paoTempSegs[0]) && (!paoTempSegs[1]))
				{
					poTempSeg = new TempSegment();
					poTempSeg->m_nMinY = i / nWidth;
					poTempSeg->m_nSubSegID = panSubSegPixels[i];
					poTempSeg->m_nSuperSegID = panSuperSegPixels[i];
					setTempSegs.insert(poTempSeg);
				}
				else if (paoTempSegs[0] && paoTempSegs[1] && (paoTempSegs[0] != paoTempSegs[1]))
				{
					poTempSeg = paoTempSegs[0];
					poTempSeg->MergeSegment(paoTempSegs[1]);
					setTempSegs.erase(paoTempSegs[1]);
					delete(paoTempSegs[1]);
				}
				else
				{
					poTempSeg = paoTempSegs[0] ? paoTempSegs[0] : paoTempSegs[1];
				}
				poTempSeg->m_setPixels.insert(i);
				poTempSeg->UpdateNeighbours(panNeighboursIndex, panSubSegPixels, panSuperSegPixels, m_mapSubSegments);


			}

			if ((i + 1) % nWidth == 0)
			{
				set<TempSegment*> setTempSegsToDelete;
				for (auto poTempSeg : setTempSegs)
				{
					if ((i + 1 == nTotalPixels) ||
						((i / nWidth) - poTempSeg->m_nMinY > SegmentMixerCatalogue::nMaxTempSegmentPixelSize))
					{
						nIncrementedMaxID = poTempSeg->CloseTempSegment(panOutputPixels, nIncrementedMaxID);
						setTempSegsToDelete.insert(poTempSeg);
					}
				}

				for (auto poTempSeg : setTempSegsToDelete)
				{
					setTempSegs.erase(poTempSeg);
					delete(poTempSeg);
				}
			}
		}


		bool bResult = (CPLE_None == poOutputDS->RasterIO(
			GF_Write, 0, 0, nWidth, nHeight, panOutputPixels, nWidth, nHeight, GDT_UInt32, 1, 0, 0, 0, 0));

		delete[]panOutputPixels;

		return bResult;
	};

	~SegmentMixerCatalogue()
	{
		for (auto el : m_mapSubSegments)
			delete(el.second);
		//for (auto el : m_mapSuperSegments)
		//	delete(el.second);

	}

protected:
	map<unsigned int, SubSegmentMeta*> m_mapSubSegments;
	//map<unsigned int, set<unsigned int>> m_mapSuperSegments;

	static const unsigned int nMaxTempSegmentPixelSize = 100;
};

/*
class SegmentMixerCatalogue : public SegmentCatalogue
{
public:
	static unsigned int* ReadDatasetIntoMem(string strSegmentsFile)
	{
		GDALDataset* poSegmentedDS = (GDALDataset*)GDALOpen(strSegmentsFile.c_str(), GA_ReadOnly);
		if (!poSegmentedDS) return false;

		int nXSize = poSegmentedDS->GetRasterXSize();
		int nYSize = poSegmentedDS->GetRasterYSize();

		unsigned int* panOutput = new unsigned int[nXSize*nYSize];

		poSegmentedDS->RasterIO(GF_Read, 0, 0, nXSize, nYSize, panOutput, nXSize, nYSize, GDT_UInt32, 1, 0, 0, 0, 0);
		GDALClose(poSegmentedDS);

		return panOutput;
	}


	//ToDO

public:
	bool InitMixedSegmentation(string strCropSegments, string strCropSubSegments)
	{
		//loop through all pixels of both files
		//newID = ordinal number of super segment * 100 000 + subsegment_id

		GDALDataset* poSegmentedDS = (GDALDataset*)GDALOpen(strCropSegments.c_str(), GA_ReadOnly);
		if (!poSegmentedDS) return false;
		int nXSize = poSegmentedDS->GetRasterXSize();
		int nYSize = poSegmentedDS->GetRasterYSize();
		unsigned int nTotal = nXSize * nYSize;
		GDALClose(poSegmentedDS);

		unsigned int* panSegmentPixels = ReadDatasetIntoMem(strCropSegments);
		unsigned int* panSubSegmentPixels = ReadDatasetIntoMem(strCropSubSegments);

		m_mapUniqueSegs.clear();
		m_mapUniqueSubSegs.clear();
		m_mapUniqueSegs[0] = 0;
		m_mapUniqueSubSegs[0] = 0;
		for (unsigned int i = 0; i < nTotal; i++)
		{
			if (m_mapUniqueSegs.find(panSegmentPixels[i]) == m_mapUniqueSegs.end())
			{
				m_mapUniqueSegs[panSegmentPixels[i]] = m_mapUniqueSegs.size();
			}
			if (m_mapUniqueSubSegs.find(panSubSegmentPixels[i]) == m_mapUniqueSubSegs.end())
			{
				m_mapUniqueSubSegs[panSubSegmentPixels[i]] = m_mapUniqueSubSegs.size();
			}
		}


		unsigned int panNeighbours[4];
		unsigned int nInd;
		unsigned int nMixedID;
		unsigned nCropStatus;
		for (unsigned int nY = 0; nY < nYSize; nY++)
		{
			for (unsigned int nX = 0; nX < nXSize; nX++)
			{

				nInd = nY * nXSize + nX;
				panNeighbours[0] = nY==0 ? 0
					: CalcMixedID(panSegmentPixels[nInd- nXSize], panSubSegmentPixels[nInd- nXSize]);
				panNeighbours[2] = nY==nYSize-1 ? 0
					: CalcMixedID(panSegmentPixels[nInd + nXSize], panSubSegmentPixels[nInd + nXSize]);
				panNeighbours[1] = nX==0 ? 0
					: CalcMixedID(panSegmentPixels[nInd - 1], panSubSegmentPixels[nInd -1]);
				panNeighbours[3] = nX == nXSize-1 ? 0
					: CalcMixedID(panSegmentPixels[nInd + 1], panSubSegmentPixels[nInd + 1]);

				nMixedID = CalcMixedID(panSegmentPixels[nInd], panSubSegmentPixels[nInd]);
				nCropStatus = panSubSegmentPixels[nInd] == 0 ? 0 : 1;

				CollectPrimaryInfoInPixel(nMixedID, nX, nY, nCropStatus, panNeighbours);
			}
		}


		delete[]panSegmentPixels;
		delete[]panSubSegmentPixels;
		return true;
	};

	//bool RunBatchProcessingOfMixedSegmentation();
	//bool AnalyzeGroupOfMixedSegments();

	bool SaveMixedSegmentation(string strOutputFile, string strCropSegments, string strCropSubSegments)
	{
		GDALDataset* poSegmentedDS = (GDALDataset*)GDALOpen(strCropSegments.c_str(), GA_ReadOnly);
		if (!poSegmentedDS) return false;
		int nXSize = poSegmentedDS->GetRasterXSize();
		int nYSize = poSegmentedDS->GetRasterYSize();
		double dblGeotransform[6];
		poSegmentedDS->GetGeoTransform(dblGeotransform);
		const char* strProjRef = poSegmentedDS->GetProjectionRef();
		//OGRSpatialReference oSRS;
		//oSRS.SetFromUserInput(strProjRef);

		unsigned int* panSegmentPixels = ReadDatasetIntoMem(strCropSegments);
		unsigned int* panSubSegmentPixels = ReadDatasetIntoMem(strCropSubSegments);
		unsigned int* panOutputPixels = new unsigned int[nXSize*nYSize];

		unsigned int nTotal = nXSize * nYSize;
		std::map<unsigned int, unsigned int> mapUniqueVals;
		mapUniqueVals[0] = 0;
		unsigned int nMixedID;
		for (unsigned int i = 0; i < nTotal; i++)
		{
			nMixedID = CalcMixedID(panSegmentPixels[i], panSubSegmentPixels[i]);
			if (mapUniqueVals.find(nMixedID) == mapUniqueVals.end())
				mapUniqueVals[nMixedID] = mapUniqueVals.size();
			panOutputPixels[i] = mapUniqueVals[nMixedID];
		}

		char **papszOptions = NULL;
		papszOptions = CSLSetNameValue(papszOptions, "COMPRESS", "LZW");

		GDALDataset* poOutputDS = (GDALDataset*)GDALCreate(
			GDALGetDriverByName("GTiff"),
			strOutputFile.c_str(),
			nXSize,
			nYSize,
			1, GDT_UInt32, papszOptions);

		poOutputDS->SetProjection(strProjRef);
		poOutputDS->SetGeoTransform(dblGeotransform);
		poOutputDS->GetRasterBand(1)->SetNoDataValue(0);

		poOutputDS->RasterIO(GF_Write, 0, 0, nXSize, nYSize, panOutputPixels, nXSize, nYSize, GDT_UInt32, 1, 0, 0, 0, 0);

		CSLDestroy(papszOptions);
		delete[]panSegmentPixels;
		delete[]panSubSegmentPixels;
		delete[]panOutputPixels;
		GDALClose(poOutputDS);
		GDALClose(poSegmentedDS);
		return true;

	};

	bool ResampleIDs()
	{
		return true;
	};

protected:
	unsigned int CalcMixedID(unsigned int nSegID, unsigned int nSubSegID)
	{
		return 100000 * m_mapUniqueSegs[nSegID] + m_mapUniqueSubSegs[nSubSegID];
	}

protected:
	std::map<unsigned int, unsigned int> m_mapUniqueSegs;
	std::map<unsigned int, unsigned int> m_mapUniqueSubSegs;
};
*/

const list<MPLOptionDescriptor> listDescriptors = {
	{"-is",			0, 0, 1, "input segmentation file" },
	{"-iss",		0, 0, 1, "input subsegmentation file" },
	{ "-o",			0, 0, 1, "output file" }
};


const list<string> listUsageExamples = {
  "adjust_segs -is segs.tif -iss sub_segs.tif -o adjust_subseg.tif",
};






#ifdef WIN32
int main(const int nArgs, const char* argv[])
{

	std::vector<string> vecArgs;
	for (int i = 0; i < nArgs; i++)
	{
		//vecArgs.push_back(argv[i]);
		vecArgs.push_back(MPLString::ReplaceAll(argv[i], "\\", "/"));
	}

	if (!MPLGDALDelayLoader::Load(MPLFileSys::RemoveExtension(vecArgs[0]) + ".config"))
	{
		cout << "ERROR: can't load GDAL" << endl;
		return 1;
	}

	cout << endl;

#else
int main(int nArgs, char* argv[])
{
	std::vector<string> vecArgs;
	for (int i = 0; i < nArgs; i++)
		vecArgs.push_back(argv[i]);
	GDALAllRegister();
	OGRRegisterAll();
	CPLSetConfigOption("OGR_ENABLE_PARTIAL_REPROJECTION", "YES");
#endif

	if (nArgs == 1)
	{
		MPLOptionParser::PrintUsage(listDescriptors, listUsageExamples);
		return 0;
	}

	MPLOptionParser oOptionParser;
	if (!oOptionParser.Init(listDescriptors, vecArgs))
	{
		cout << "ERROR: input cmd line is not valid" << endl;
		return 1;
	}

	SegmentMixerCatalogue oCatalogue;

	unsigned int nWidth, nHeight;
	unsigned int* panSubSegPixels = (unsigned int*)ReadAllPixelsFromSingleBandRaster(oOptionParser.GetOptionValue("-iss"),
		nWidth, nHeight);
	unsigned int* panSuperSegPixels = (unsigned int*)ReadAllPixelsFromSingleBandRaster(oOptionParser.GetOptionValue("-is"),
		nWidth, nHeight);
	//debug
	//nHeight = 500;
	//end-debug
	oCatalogue.InitFromArrays(panSubSegPixels, panSuperSegPixels, nWidth, nHeight);
	//oCatalogue.CalcSubstitutes();

	GDALDataset* poOutputDS = CreateGDALDatasetForWritingUsingPattern(oOptionParser.GetOptionValue("-o"),
		oOptionParser.GetOptionValue("-iss"));
	oCatalogue.SaveAdjustedSubsegments(panSubSegPixels, panSuperSegPixels, nWidth, nHeight, poOutputDS);

	//debug
	std::cout << nNoSubstituteCount << std::endl;
	std::cout << nSubstituteCount << std::endl;
	//end-deubug
	GDALClose(poOutputDS);



	/*
	oCatalogue.InitMixedSegmentation(oOptionParser.GetOptionValue("-is"), oOptionParser.GetOptionValue("-iss"));
	std::cout << oCatalogue.GetSize() << endl;
	oCatalogue.SaveMixedSegmentation(oOptionParser.GetOptionValue("-o"),
		oOptionParser.GetOptionValue("-is"),
		oOptionParser.GetOptionValue("-iss"));
	*/

	return 0;
}