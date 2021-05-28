// extract_fields.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
// extract_fields.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include "external.h"
#include "merge_segs.h"


const list<MPLOptionDescriptor> listDescriptors = {
	{"-ic",			0, 0, 1, "input crop pixel mask" },
	{"-is",			0, 0, 1, "input segments file" },
	{"-ie",			0, 0, 1, "input edges file" },
	{ "-o",			0, 0, 1, "output file" },
	{"-mina",		0, 0, 0, "min seg. area in pixels (default 0)"},
	{"-merge",		1, 0, 0, "try merge segments "},
	{"-maxh",		0, 0, 0, "max border height threshold (default 0.35) " },
	{"-minh",		0, 0, 0, "max border height threshold (default 0.20) " },
	{"-del_small",	1, 0, 0, "delete small not merged crop segments"},
	{"-del_marg",	1, 0, 0, "delete marginal segments"},
	{"-add",		0, 0, 0, "add value to separate segments in different tiles"}
	//	{"-incl",	0, 0, 0, "inclusiveness threshold  (default 0.66)" }
		//{"-drop_no_crop", 1, 0, 0, "delete all none crop segments"}
};


const list<string> listUsageExamples = {
  "merge_segs -ic crop_mask_map_39UVB.tif -is segmentation-merged-39UVB.tif -o result.tif -m 10",
};


const int DEFAULT_MIN_AREA = 0;
const float DEFAULT_MAXHEIGHT = 0.2;
const float DEFAULT_MINHEIGHT = 0.35;
const float DEFAULT_INCLUSIVENESS = 0.66;




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

	std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();


	cout << "1. Gathering segment metadata ... ";
	SegmentCatalogue oCatalogue;
	if (!oCatalogue.InitFromFiles(oOptionParser.GetOptionValue("-is"),
		oOptionParser.GetOptionValue("-ic")))
	{
		std::cout << "ERROR: SegmentCatalogue::InitFromFiles fail" << std::endl;
		return 2;
	}

	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
	std::cout << "done" << std::endl;
	std::cout << "Segments count before processing: " << oCatalogue.GetSize() << endl;
	std::cout << "Duration = "
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;

	int nMinArea = oOptionParser.GetOptionValue("-mina") == "" ? DEFAULT_MIN_AREA
		: std::atoi(oOptionParser.GetOptionValue("-mina").c_str());

	unsigned int nLastNotRemovedID = 0;

	//debug
	//SegmentMeta* poSeg = oCatalogue.GetSegmentRef(674073);
	//std::cout << poSeg->nArea << endl;
	//std::cout << poSeg->nBorderLength << endl;
	//std::cout << poSeg->nCropPixelsCount << endl;
	//end-debug

	//Primary segment filtering (first segment postprocessing procedure) -
	//using simple condition on crop pixels/area values 

	begin = std::chrono::steady_clock::now();
	std::cout << endl;
	std::cout << "2. Deleting none crop segments ... ";
	oCatalogue.DropAllNoCrop();
	std::cout << "done" << std::endl;
	std::cout << "Crop segments: " << oCatalogue.GetSize() << endl;
	end = std::chrono::steady_clock::now();
	std::cout << "Duration = "
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;


	if (oOptionParser.GetOptionValue("-merge") != "" ||
		oOptionParser.GetOptionValue("-mina") != "")
	{
		std::cout << endl;
		cout << "3. Merging and filtering crop segments ... ";
		begin = std::chrono::steady_clock::now();

		if (oOptionParser.GetOptionValue("-merge") != "")
		{
			if (oOptionParser.GetOptionValue("-ie") == "")
			{
				std::cout << "ERROR: -ie input isn't specified" << endl;
				return 0;
			}

			oCatalogue.InitFromFiles(oOptionParser.GetOptionValue("-is"),
				oOptionParser.GetOptionValue("-ic"),
				oOptionParser.GetOptionValue("-ie"));

			float dblMaxHeight = oOptionParser.GetOptionValue("-maxh") == "" ? DEFAULT_MAXHEIGHT
				: std::atof(oOptionParser.GetOptionValue("-maxh").c_str());

			float dblMinHeight = oOptionParser.GetOptionValue("-minh") == "" ? DEFAULT_MINHEIGHT
				: std::atof(oOptionParser.GetOptionValue("-minh").c_str());

			//float dblInclussiveness = oOptionParser.GetOptionValue("-incl") == "" ? DEFAULT_INCLUSIVENESS
			//					: std::atof(oOptionParser.GetOptionValue("-incl").c_str());

			oCatalogue.RunBatchMerging(nMinArea, dblMinHeight, dblMaxHeight);
			std::cout << oCatalogue.GetSize() << " ";
			oCatalogue.MergeIfInsideCropGroup(nMinArea);
			std::cout << oCatalogue.GetSize() << " ";
			oCatalogue.RunBatchMerging(nMinArea, dblMinHeight, dblMaxHeight);
			std::cout << oCatalogue.GetSize() << " ";
		}

		if (oOptionParser.GetOptionValue("-del_small") != "")
		{
			oCatalogue.FilterNotMergedSmallSegments(nMinArea);
		}

		std::cout << "done " << std::endl;
		std::cout << "Final count after filtering segments: " << oCatalogue.GetSize() << endl;
		end = std::chrono::steady_clock::now();
		std::cout << "Duration = "
			<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;

	}
	else cout << "3. Missing filtering step" << endl;


	std::cout << endl;
	cout << "4. Saving output ... ";
	begin = std::chrono::steady_clock::now();
	unsigned int nAddVal = oOptionParser.GetOptionValue("-add") == "" ? 0
		: std::atoi(oOptionParser.GetOptionValue("-add").c_str());
	oCatalogue.SaveMergedSegments(oOptionParser.GetOptionValue("-o"),
		oOptionParser.GetOptionValue("-is"),
		oOptionParser.GetOptionValue("-del_marg") != "",
		nAddVal);
	std::cout << "done " << std::endl;
	end = std::chrono::steady_clock::now();
	std::cout << "Duration = "
		<< std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "s" << std::endl;
	return 0;
}


