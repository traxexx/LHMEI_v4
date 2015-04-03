#include <string> 
#include <dirent.h>
#include <iostream>
#include <utility>

#include "ComputeLHMEI.h"
#include "Utilities.h" // file check , done file generate etc
#include "Preprocess.h" // generate ref-bam & disc_bam
#include "ReadMap.h"
#include "bamUtility.h" // getAvrReadLen;
#include "Globals.h" // DEBUG_MODE, WIN, STEP

#include "RefStats.h"

using std::string;
using std::cout;
using std::cerr;
using std::endl;

void ComputeLHMEI (Options * ptrMainOptions)
{
	int win_len = stoi(ptrMainOptions->ArgMap["Win"]);
	int step_len = stoi(ptrMainOptions->ArgMap["Step"]);
	WIN = stoi(ptrMainOptions->ArgMap["Win"]);
	STEP = stoi(ptrMainOptions->ArgMap["Step"]);
	string REF_CHR = ptrMainOptions->ArgMap["CtrlChr"];
	if (ptrMainOptions->OptMap["debug"])
		DEBUG_MODE = 1;

	string work_dir = ptrMainOptions->ArgMap["WorkDir"];
	if (work_dir[ work_dir.size() - 1 ] != '/')
		work_dir += '/';
	string bam_dir = work_dir + "bam_tmps/";
	string ctrl_dir = work_dir + "ctrl_tmps/";
	string pre_dir = work_dir + "preprocess/";

// prepare directories	
	string cmd = "mkdir -p " + bam_dir + " " + ctrl_dir + " " + pre_dir;
	ExecuteCmd(cmd);	
	
/*** pre-process: do not focus on any chr ***/
	string outSam = pre_dir + REF_CHR + ".sam";
	string discNsortBam = pre_dir + "disc-nsort";
	string discSam = pre_dir + "disc.sam";
	int avr_read_len = GetAvrReadLenFromBam( ptrMainOptions->ArgMap["Bam"].c_str() );
	int avr_ins_size = GetAvrInsSizeFromBam( ptrMainOptions->ArgMap["Bam"].c_str() );
	if (!ExistDoneFile( pre_dir, "PreProcess" )) {
		PreProcessBam( ptrMainOptions->ArgMap["Bam"].c_str(), outSam.c_str(), REF_CHR.c_str(),
			discSam.c_str(), ptrMainOptions->ArgMap["MEcoord"].c_str(),
			avr_read_len, avr_ins_size);
		string disc_nsort_cmd = string("samtools sort -n ") + discSam + " " + discNsortBam;
		ExecuteCmd(disc_nsort_cmd);
		GenerateDoneFile( pre_dir, "PreProcess" );
	}
	discNsortBam += ".bam";

/*** Generate level list if >= 0 ***/
  // generate bam index of discSam first
	LEVEL = stoi( ptrMainOptions->ArgMap["Simplify"] );
	string level_file_prefix = pre_dir + "level-list";
	if (LEVEL >= 0) {
		if ( !ExistDoneFile( pre_dir, "LevelList" ) ) {
			GenerateLevelListFromDiscBam( discNsortBam.c_str(), level_file_prefix.c_str(), step_len, win_len );
			GenerateDoneFile( pre_dir, "LevelList" );	
		}
	}

/*** re-map ***/
	string ctrl_bam = ctrl_dir + REF_CHR + "-remap-sort-recal.bam";
	if ( !ExistDoneFile(ctrl_dir, "Remap") ) {
	  // nsort
		string cmd = string("samtools sort -n ") + outSam + " " + ctrl_dir + REF_CHR + "-nsort";
		ExecuteCmd(cmd);
	  // fastq
	  	cmd = string("bam bam2FastQ --in ") + ctrl_dir + REF_CHR + "-nsort.bam --readName --outBase " + ctrl_dir + REF_CHR;
	  	ExecuteCmd(cmd);
	  // re-map
	  	string slice_ref = ptrMainOptions->ArgMap["ProgramDir"] + "/refs/slice-chr20-hs37d5.fa ";
	  	string fastq_prefix = ctrl_dir + REF_CHR;
	  	string remapSam = ctrl_dir + "align-pe.sam";
	  	cmd = GetRemapCmd(ptrMainOptions->ArgMap["Mapper"], fastq_prefix, slice_ref, remapSam);
	  	ExecuteCmd(cmd);
	  // sort by coord
	  	cmd = string("samtools sort ") + remapSam + " " + ctrl_dir + REF_CHR + "-remap-sort";
	  	ExecuteCmd(cmd);
	  // dedup
	  	cmd = string("bam dedup --recab --in ") + ctrl_dir + REF_CHR + "-remap-sort.bam --out ";
	  	cmd += ctrl_bam + " --force --refFile " + slice_ref + "--storeQualTag OQ --maxBaseQual 40";
	  	ExecuteCmd(cmd);
	  // generate .bai
	  	cmd = string("samtools index ") + ctrl_bam;
	  	ExecuteCmd(cmd);
	  // generate done
	  	GenerateDoneFile( ctrl_dir, "Remap" );
	}


	string focus_chr_str = ptrMainOptions->ArgMap["Chr"];
/*** Read Map (target + ctrl )****/
	string bam_done_flag = string("BamMap-");
	bool bam_counted;
	if ( focus_chr_str.compare("-1")  == 0 ) {
		bam_done_flag += string("all");
		bam_counted = ExistDoneFile( bam_dir, bam_done_flag.c_str() );
	}
	else {
		string all_done_flag = bam_done_flag; // if all done, then consider single-chr is done
		all_done_flag += string("all");
		bam_done_flag += focus_chr_str;
		bam_counted = ( ExistDoneFile( bam_dir, bam_done_flag.c_str() ) | ExistDoneFile( bam_dir, all_done_flag.c_str() ) );
	}
	bool ctrl_counted = ExistDoneFile( ctrl_dir, "CtrlBamMap" );
	if ( !bam_counted || !ctrl_counted ) {
		ReadMap * Rmap = new ReadMap( win_len, step_len, avr_read_len, avr_ins_size,
			REF_CHR.c_str(), ptrMainOptions->ArgMap["MElist"].c_str() );
		const char * focus_chr = (focus_chr_str.compare("-1") == 0) ? "" : focus_chr_str.c_str();
		if ( !bam_counted ) {
			std::cout << "Setting read types from raw bam..." << std::endl;
			Rmap->SetMapFromBam(ptrMainOptions->ArgMap["Bam"].c_str(), discNsortBam.c_str(), bam_dir.c_str(),
				focus_chr, level_file_prefix.c_str());
			GenerateDoneFile( bam_dir, bam_done_flag.c_str() );
		}
		if ( !ctrl_counted ) {
			std::cout << "Setting read types from ctrl bam..." << std::endl;
			string ctrl_coord = ptrMainOptions->ArgMap["MEcoord"] + ".liftOver";
			Rmap->SetMapFromCtrlBam( ctrl_bam.c_str(), ctrl_dir.c_str(), ctrl_coord.c_str() );
			GenerateDoneFile( ctrl_dir, "CtrlBamMap" );
		}
		delete Rmap;
	}
	
/*** Calculate LH ****/
// here do not set done file
	AllHetIndex allHetIndex;
	SetAllHetIndex( ptrMainOptions->ArgMap["HetIndex"].c_str(), allHetIndex );
	string sample_name = ptrMainOptions->ArgMap["-Sample"];
// do by mei type
	for( int mei_type = 0; mei_type <= 2; mei_type++ ) {
		cout << "Discovering mei-type: " << mei_type << " ..." << endl;
		string lh_flag = string("LH-") + focus_chr_str + "." + std::to_string(mei_type);
		if ( !ExistDoneFile( bam_dir, lh_flag.c_str() ) ) {	
			string ctrl_proper_prefix = ctrl_dir + "proper-" + REF_CHR;
			string ctrl_disc_prefix = ctrl_dir + "disc-" + REF_CHR;
		// construct ref stats
			RefStats* rStats = new RefStats( ctrl_proper_prefix, ctrl_disc_prefix, mei_type, allHetIndex );
			string outRecord = work_dir + "refLH." + std::to_string(mei_type) + ".report";
		// ref LH
			rStats->PrintCtrlGLasRecord( outRecord, 1 );
			rStats->MarkRefLHasDone();
			rStats->ReAdjustSelfsWithLiftOver();
		// data LH: loop-through bam dir ---> re-organize --> 
			OriginalStats* dataOsPtr = new OriginalStats( mei_type, sample_name );
		  //add
			if ( focus_chr_str.compare("-1") != 0 ) { // single chr
				string data_proper_name = bam_dir + "proper-" + focus_chr_str;
				string data_disc_name = bam_dir + "disc-" + focus_chr_str;
				dataOsPtr->Add( focus_chr_str, data_proper_name, data_disc_name );
			}
			else { // loop through bam dir
				struct dirent *pDirent;
				DIR *pDir;
				pDir = opendir( bam_dir.c_str() );
				if ( pDir == NULL ) {
					cerr << "ERROR: Cannot open " << bam_dir << endl;
					exit(1);
				}
				while ((pDirent = readdir(pDir)) != NULL) {
					string data_proper_name = pDirent->d_name;
					if ( data_proper_name[ data_proper_name.size() - 1 ] != '3' ) // only read 3 since disc do not have 3
						continue;
					data_proper_name.substr(0, data_proper_name.size() - 2);
					int slash_loc = 0;
					for( int i = data_proper_name.size() - 1; i >= 0; i-- ) {
						if ( data_proper_name[i] == '/' ) {
							slash_loc = i;
							break;
						}
					}
					if ( slash_loc == 0 ) {
						cerr << "ERROR: no dir name..." << endl;
						exit(1);
					}
					int first_dot = 0;
					int second_dot = 0;
					for( unsigned int i = slash_loc + 1; i < data_proper_name.size(); i++) {
						if ( data_proper_name[i] == '.' ) {
							if (first_dot > 0) {
								second_dot = i;
								break;
							}
							else
								first_dot = i;
						}
					}
					if ( second_dot == 0 ) {
						cerr << "ERROR: no dot..." << endl;
						exit(1);
					}
					
					string current_chr = data_proper_name.substr( first_dot + 1, second_dot - first_dot );				
					data_proper_name = data_proper_name.substr( 0, data_proper_name.size() - 2 );
					string data_disc_name = data_proper_name;
					data_disc_name.replace( slash_loc + 1, 6, string("disc") );
					dataOsPtr->Add( current_chr, data_proper_name, data_disc_name );
				}
				closedir( pDir );
			}
		// re-organize & clear under level
			dataOsPtr->ReOrganize();
			dataOsPtr->ClearUnderLevelMergeCells();
		// calculate LH
			for( MergeCellPtr merge_it = dataOsPtr->MergeData.begin(); merge_it != dataOsPtr->MergeData.end(); merge_it++ ) {
			  	rStats->SetRecordGL( merge_it );
			}
			rStats->AdjustUsedLoci( dataOsPtr );
			string vcf_name = work_dir + "Hits";
			if ( focus_chr_str.compare("-1") != 0 )
				vcf_name += "-" + focus_chr_str;	
			 vcf_name += "." + std::to_string(mei_type) + ".vcf";
			dataOsPtr->PrintGLasVcf( vcf_name );
		}
	}
// delete intermediates
}




















