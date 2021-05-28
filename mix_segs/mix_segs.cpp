// mix_segs.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../merge_segs/merge_segs.h"



struct SegmentMetaMultiLevel : public SegmentMeta
{

	SegmentMetaMultiLevel* m_poParentSeg;
	map<unsigned int, SegmentMetaMultiLevel*> m_mapChildSegs;
	bool m_bIsSplitted;

	SegmentMetaMultiLevel()
	{
		m_poParentSeg = 0;
		m_bIsSplitted = false;
	};

	~SegmentMetaMultiLevel()
	{
		//Todo may be some recursive 
	};

	bool IntiMultilevelSegmentationFromInMemData(unsigned int nWidth,
		unsigned int nHeight,
		uint16_t* panCropMaskPixels,
		vector<unsigned int*> vectSegPixels
	)
	{
		unsigned int nTotalPixels = nWidth*nHeight;
		SegmentMetaMultiLevel* poMLSeg = 0;
		unsigned int nLevels = vectSegPixels.size();
		vector<unsigned int> vectOnePixelByLevel(vectSegPixels.size(), 0);
		vector<vector<unsigned int>> vectNeighboursByLevel(vectSegPixels.size(), {0,0,0,0});

		for (unsigned int i = 0; i < nTotalPixels; i++)
		{
			for (unsigned int l = 0; l < nLevels; l++)
			{
				vectOnePixelByLevel[l] = vectSegPixels[l][i];
				vectNeighboursByLevel[l][0] = i < nWidth ? 0 : vectSegPixels[l][i-nWidth];
				vectNeighboursByLevel[l][1] = i % nWidth == 0 ? 0 : vectSegPixels[l][i - 1];
				vectNeighboursByLevel[l][2] = i + nWidth >= nTotalPixels ? 0 : vectSegPixels[l][i + nWidth];
				vectNeighboursByLevel[l][3] = (i + 1) % nWidth == 0 ? 0 : vectSegPixels[l][i + 1];
			}
			
			SegmentMetaMultiLevel::ProcessPixelRecursive(this,
				panCropMaskPixels[i],
				0,
				nLevels,
				vectOnePixelByLevel,
				vectNeighboursByLevel);
		}

		return true;
	}

	//static unsigned int CalcOutputPixelRecursive(SegmentMetaMultiLevel* poMLSeg, uint16_t )
	//parse input rasters
//store all segs into some tree type data structure
//loop through level 1 segs and calc output values

//general scheme:
// we should analyze lev1 taking into account all subsegs including lev4
// if lev 1 isn't satisfied than we take take first sub level 2 or 3 where there are several subsegments and than repeat analysis

//rules:
//should split anyway if area is large than 400 pixels
//should split crop seg if:
//		sum of none-crop segs (lev1) is greater than 2 ha or larger than 5 percent of area
//should split none-crop seg if
//		there is at least one crop subseg greater than 2 ha
//		or sum of area of crop segs is greater than
//after those steps only topological cases are left

//ToDo
//1. Update CalcNoneCropSubSegsArea and CalcCropSubSegsArea
//   - should take into account large subsegs if there are exists than split
//
//1b. Split large segs ??? may be not a good idea
//
//2. Solve problem of large crop segs that contains several neighbouring subsegs
//Algorithm:
// if there are not more than 3 subsegs that cover more than 90-95 of area than split
//
//
//3. Merging algorithm problem
//
//4. Fix nID for output
//
//

	void CalcNoneCropSubSegsArea(unsigned int &nNoneCropSubSegsArea)
	{
		for (auto poMLSubSeg : m_mapChildSegs)
		{
			if (!poMLSubSeg.second->IsRoughlyCrop())
				nNoneCropSubSegsArea += poMLSubSeg.second->nArea;
			else
				poMLSubSeg.second->CalcNoneCropSubSegsArea(nNoneCropSubSegsArea);
		}
	}

	void CalcCropSubSegsArea(unsigned int &nCropSubSegsArea)
	{
		for (auto poMLSubSeg : m_mapChildSegs)
		{
			if (poMLSubSeg.second->IsRoughlyCrop())
				nCropSubSegsArea += poMLSubSeg.second->nArea;
			else
				poMLSubSeg.second->CalcCropSubSegsArea(nCropSubSegsArea);
		}
	}

	bool SplitIfFewLargeSubSegs(unsigned int nMinArea, unsigned int nMaxLargeSubSegs, double dfAreaCoeff)
	{
		//bypass all subsegs and select three largest
		//
		//
		for (auto poSubSeg : m_mapChildSegs)
		{
			if (poSubSeg.second->nArea < nMinArea) continue;
			if (poSubSeg.second->m_bIsSplitted)
			{
				SegmentMetaMultiLevel::SplitIfFewLargeSubSegsRecursive(poSubSeg.second,
					nMinArea, nMaxLargeSubSegs, dfAreaCoeff);
			}
			else if (poSubSeg.second->IsRoughlyCrop())
			{
				if (poSubSeg.second->m_mapChildSegs.size() == 1) continue; //Not accurate, to be fixed further
				else SegmentMetaMultiLevel::SplitIfFewLargeSubSegsRecursive(poSubSeg.second,
					nMinArea, nMaxLargeSubSegs, dfAreaCoeff);
			}
		}

		


		return true;
	}


	bool SplitIfOppositeTypeSubSegs(unsigned int nMinArea, unsigned int nMaxArea)
	{
		unsigned int nCropSubSegsArea;
		unsigned int nNoneCropSubSegsArea;

		for (auto poSubSeg : m_mapChildSegs)
		{
			nNoneCropSubSegsArea = 0;
			nCropSubSegsArea = 0;
			if (poSubSeg.second->nArea < nMinArea) continue;
			
			if (poSubSeg.second->IsRoughlyCrop())
			{
				if (poSubSeg.second->nArea > nMaxArea)
				{
					poSubSeg.second->m_bIsSplitted = true;
				}
				else
				{
					poSubSeg.second->CalcNoneCropSubSegsArea(nNoneCropSubSegsArea);
					if (20 * nNoneCropSubSegsArea > poSubSeg.second->nArea)
					{
						poSubSeg.second->m_bIsSplitted = true;
					}

				}
			}
			else
			{
				poSubSeg.second->CalcCropSubSegsArea(nCropSubSegsArea);
				if (20 * nCropSubSegsArea > poSubSeg.second->nArea)
				{
					poSubSeg.second->m_bIsSplitted = true;
				}
			}

			if (poSubSeg.second->m_bIsSplitted && (poSubSeg.second->m_mapChildSegs.size()>1))
				poSubSeg.second->SplitIfOppositeTypeSubSegs(nMinArea, nMaxArea);
			
		}

		return true;
	}

	unsigned int GetLastLevelID()
	{
		if (m_mapChildSegs.size() == 0) return nID;
		else return (*m_mapChildSegs.begin()).second->GetLastLevelID();
	}

	bool SaveOutput(GDALDataset* poOutputDS,
		unsigned int nWidth,
		unsigned int nHeight,
		vector<unsigned int*> vectSegPixels
	)
	{
		unsigned int nTotalPixels = nWidth * nHeight;

		SegmentMetaMultiLevel* poMLSeg = 0;
		unsigned int nLevels = vectSegPixels.size();
		unsigned int* panOutputPixels = new unsigned int[nTotalPixels];

		for (int i = 0; i < nTotalPixels; i++)
		{
			poMLSeg = this;
			for (unsigned int l = 0; l < nLevels; l++)
			{
				poMLSeg = poMLSeg->m_mapChildSegs[vectSegPixels[l][i]];
				if (!poMLSeg->m_bIsSplitted)
				{

					panOutputPixels[i] = poMLSeg->GetLastLevelID();// poMLSeg->nID;
					break;
				}
			}
		}

		poOutputDS->RasterIO(GF_Write, 0, 0, nWidth, nHeight, panOutputPixels, nWidth, nHeight, GDT_UInt32, 1, 0, 0, 0, 0, 0);
		

		delete[]panOutputPixels;
		return true;
	}

	static void CalcStatRecursive(SegmentMetaMultiLevel* poMLSeg,
		const unsigned int &nMaxLevel,
		unsigned int nLevel,
		unsigned int &nMaxOneLevelSplit,
		unsigned int &nSegIDWithMaxSplit)
	{

		if (poMLSeg->m_mapChildSegs.size() > nMaxOneLevelSplit)
		{
			nMaxOneLevelSplit = poMLSeg->m_mapChildSegs.size();
			nSegIDWithMaxSplit = poMLSeg->nID;
		}
		if (nLevel + 2 == nMaxLevel) return;

		for (auto poMLChildSeg : poMLSeg->m_mapChildSegs)
		{
			CalcStatRecursive(poMLChildSeg.second, nMaxLevel, nLevel + 1, nMaxOneLevelSplit, nSegIDWithMaxSplit);
		}

	}

	static bool CalclStat(SegmentMetaMultiLevel* poParentMLSeg, unsigned int nMaxLevel)
	{
		//num of Lev0 segs without split
		//segs with max number subsegs total and one lev split
		unsigned int nNoSplitSegs = 0;
		unsigned int nSegIDWithMaxSplit = 0;
		unsigned int nMaxOneLevelSplit = 0;
		SegmentMetaMultiLevel* poMLSeg = 0;

		for (auto poLev0MLSeg : poParentMLSeg->m_mapChildSegs)
		{
			unsigned int l;
			poMLSeg = poLev0MLSeg.second;
			for (l = 0; l < nMaxLevel-1; l++)
			{
				if (poMLSeg->m_mapChildSegs.size() != 1) break;
				else poMLSeg = (*poMLSeg->m_mapChildSegs.begin()).second;
			}
			if (l == (nMaxLevel-1)) nNoSplitSegs++;
		
			if ((poLev0MLSeg.first!=1) && (2*poLev0MLSeg.second->nCropPixelsCount > poLev0MLSeg.second->nArea)) 
				CalcStatRecursive(poLev0MLSeg.second, nMaxLevel, 0, nMaxOneLevelSplit, nSegIDWithMaxSplit);

		}

		std::cout << nNoSplitSegs << std::endl;
		std::cout << nMaxOneLevelSplit << std::endl;
		std::cout << nSegIDWithMaxSplit << std::endl;

		return true;
	}

	static bool SplitIfFewLargeSubSegsRecursive(SegmentMetaMultiLevel* poMLSeg,
		unsigned int nMinArea,
		unsigned int nMaxLargeSubSegs,
		double dfAreaCoeff)
	{
		if (poMLSeg->nArea < nMinArea) return true;
		else if (poMLSeg->m_mapChildSegs.size() == 1) return true;
		else
		{
			list<unsigned int> listSubSegsArea;
			for (auto poSubSeg : poMLSeg->m_mapChildSegs)
			{
				listSubSegsArea.push_back(poSubSeg.second->nArea);
			}

			listSubSegsArea.sort(std::greater<unsigned int>());

			int i = 0;
			unsigned int nAreaSum = 0;
			for (auto el : listSubSegsArea)
			{
				nAreaSum += el;
				i++;
				if (i == nMaxLargeSubSegs) break;
			}
			//debug
			if (poMLSeg->nID == 742601)
			{
				std::cout << poMLSeg->m_mapChildSegs.size() << std::endl;
				std::cout << poMLSeg->nArea << " " << nAreaSum << " " << (unsigned int)(dfAreaCoeff * poMLSeg->nArea) << std::endl;
			}
			//end-debug

			if (nAreaSum > (unsigned int)(dfAreaCoeff * poMLSeg->nArea))
			{
				poMLSeg->m_bIsSplitted = true;
				return true; //Not honest
				for (auto poSubSeg : poMLSeg->m_mapChildSegs)
					SegmentMetaMultiLevel::SplitIfFewLargeSubSegsRecursive(poSubSeg.second,
						nMinArea, nMaxLargeSubSegs, dfAreaCoeff);
			}
		}

		return true;
	}

	static void ProcessPixelRecursive(SegmentMetaMultiLevel* poParentMLSeg,
									uint16_t nCropStatus,
									unsigned int nLevel,
									unsigned int nMaxLevel,
									vector<unsigned int> &vectOnePixelByLevel,
									vector<vector<unsigned int>> &vectNeighboursByLevel
									)
	{
		SegmentMetaMultiLevel* poNewParentMLSeg = 
			poParentMLSeg->m_mapChildSegs.find(vectOnePixelByLevel[nLevel]) == poParentMLSeg->m_mapChildSegs.end()
				? new SegmentMetaMultiLevel() : poParentMLSeg->m_mapChildSegs[vectOnePixelByLevel[nLevel]];

		poNewParentMLSeg->nID = poNewParentMLSeg->nID == 0 ? vectOnePixelByLevel[nLevel] :poNewParentMLSeg->nID;
		poNewParentMLSeg->nArea += 1;
		poNewParentMLSeg->nCropPixelsCount += nCropStatus;

		///*
		for (int i = 0; i < 4; i++)
		{
			if (vectNeighboursByLevel[nLevel][i] && (vectNeighboursByLevel[nLevel][i] != poNewParentMLSeg->nID))
			{
				poNewParentMLSeg->nBorderLength++;
				break;
			}
		}
		//*/

		poParentMLSeg->m_mapChildSegs[vectOnePixelByLevel[nLevel]] = poNewParentMLSeg;

		if (nLevel + 1 < nMaxLevel) 
			ProcessPixelRecursive(poNewParentMLSeg, 
				nCropStatus, nLevel + 1, nMaxLevel, vectOnePixelByLevel, vectNeighboursByLevel);

		return;
	}
};


//parse input rasters
//store all segs into some tree type data structure
//loop through level 1 segs and calc output values

//general scheme:
// we should analyze lev1 taking into account all subsegs including lev4
// if lev 1 isn't satisfied than we take take first sub level 2 or 3 where there are several subsegments and than repeat analysis

//rules:
//should split anyway if area is large than 400 pixels
//should split crop seg if:
//		sum of none-crop segs (lev1) is greater than 2 ha or larger than 5 percent of area
//should split none-crop seg if
//		there is at least one crop subseg greater than 2 ha
//		or sum of area of crop segs is greater than
//after those steps only topological cases are left


const unsigned int DEF_MIN_AREA = 200;
const unsigned int DEF_MAX_AREA = 35000;
const unsigned int DEF_LARGE_SUB_SEGS = 3;
const unsigned int DEF_AREA_COEFF = 0.95;


const list<MPLOptionDescriptor> listDescriptors = {
	{"-is",			0, 1, 1, "input segmentation file" },
	{"-ic",			0, 0, 1, "input pixel crop mask file" },
	{ "-o",			0, 0, 1, "output file" },
	{ "-mina",		0, 0, 0, "min area to split (default: 200)" },
	{ "-maxa",		0, 0, 0, "min area to split (default: 35000)" }
	//{ "-ls",		0, 0, 0, "large sub. segs. (default: 3)" },
	//{ "-ac",		0, 0, 0, "area coeff. to split  (default: 0.95)" }

};


const list<string> listUsageExamples = {
  "mix_segmentations -is segs.tif -iss sub_segs.tif -o adjust_subseg.tif",
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

	


	cout << "1. Reading input pixels into mem ... ";
	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


	unsigned int nWidth, nHeight;
	uint16_t* panCropMaskPixels = (uint16_t*)ReadAllPixelsFromSingleBandRaster(
		oOptionParser.GetOptionValue("-ic"), nWidth, nHeight);
	if (!panCropMaskPixels)
	{
		std::cout << "ERROR: can't read file: " << oOptionParser.GetOptionValue("-ic") << std::endl;
		return 1;
	}

	vector<unsigned int*> vectSegPixels;
	for (auto strSegFile : oOptionParser.GetValueList("-is"))
	{
		vectSegPixels.push_back((unsigned int*)ReadAllPixelsFromSingleBandRaster(
			strSegFile,nWidth, nHeight));
		if (!vectSegPixels.back())
		{
			std::cout << "ERROR: can't read file: " << strSegFile << std::endl;
			return 2;
		}
	}
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "done in " 
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;

	SegmentMetaMultiLevel oSegMultiLevel;
	//debug
	//nHeight = 500;
	//end-debug

	cout << "2. Init. multi level segments  ... ";
	begin = std::chrono::steady_clock::now();
	oSegMultiLevel.IntiMultilevelSegmentationFromInMemData(nWidth, nHeight, panCropMaskPixels, vectSegPixels);
	end = std::chrono::steady_clock::now();
	std::cout << "done in "
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;
	
	//debug
	std::cout << oSegMultiLevel.m_mapChildSegs.size() << std::endl;
	SegmentMetaMultiLevel::CalclStat(&oSegMultiLevel, 4);

	//end-debug

	unsigned int nMinArea = oOptionParser.GetOptionValue("-mina") == "" ? DEF_MIN_AREA :
		std::atoi(oOptionParser.GetOptionValue("-mina").c_str());

	unsigned int nMaxArea = oOptionParser.GetOptionValue("-maxa") == "" ? DEF_MAX_AREA :
		std::atoi(oOptionParser.GetOptionValue("-maxa").c_str());

	//unsigned int nMaxLargeSubSegs = oOptionParser.GetOptionValue("-ls") == "" ? DEF_LARGE_SUB_SEGS :
	//	std::atoi(oOptionParser.GetOptionValue("-ls").c_str());

	//double dfAreaCoeff = oOptionParser.GetOptionValue("-ac") == "" ? DEF_AREA_COEFF :
	//	std::atof(oOptionParser.GetOptionValue("-ac").c_str());


	cout << "3. Analyze segments for splitting ... ";
	begin = std::chrono::steady_clock::now();
	oSegMultiLevel.SplitIfOppositeTypeSubSegs(nMinArea,nMaxArea);
	//oSegMultiLevel.SplitIfFewLargeSubSegs(nMinArea, nMaxLargeSubSegs, dfAreaCoeff);

	end = std::chrono::steady_clock::now();
	std::cout << "done in "
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;


	cout << "4. Saving output ... ";
	begin = std::chrono::steady_clock::now();
	GDALDataset* poOutputDS = CreateGDALDatasetForWritingUsingPattern(oOptionParser.GetOptionValue("-o"),
		(*oOptionParser.GetValueList("-is").begin()));

	if (!poOutputDS)
	{
		std::cout << "ERROR: can't create output file: " << oOptionParser.GetOptionValue("-o") << std::endl;
		return 3;
	}

	bool bResult = oSegMultiLevel.SaveOutput(poOutputDS, nWidth, nHeight, vectSegPixels);
	GDALClose(poOutputDS);
	end = std::chrono::steady_clock::now();
	std::cout << "done in "
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;
	

	return 0;
}
