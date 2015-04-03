#include <iostream>
#include <stdlib.h>
#include <utility> // exit
#include "Utilities.h"
#include <iterator>


void CheckInputFileStatus(std::ifstream & infile, const char* name)
{
	if (!infile.is_open()) {
		std::cerr << "ERROR: Can't open input file " << name << std::endl;
		exit(1);
	}
}

void CheckOutFileStatus(std::ofstream & outfile, const char* name)
{
	if (!outfile.is_open()) {
		std::cerr << "ERROR: Can't open output file " << name << std::endl;
		exit(1);
	}
}


std::string getMeiTypeString(int & mei_type)
{
	std::string MeiType;
	switch( mei_type ) {
		case 0:
			MeiType = std::string("Alu"); break;
		case 1:
			MeiType = std::string("L1"); break;
		case 2:
			MeiType = std::string("SVA"); break;
		default:
			std::cerr << "ERROR: illegal mei_type in ReadMap::PrintToVcf!" << std::endl; exit(1);
	}
	return MeiType;
}

// use as DP in vcf
int getSumOfVector( vector<int> & count_table )
{
	int sum = 0;
	for( vector<int>::iterator it = count_table.begin(); it != count_table.end(); it++ )
		sum += *it;
	return sum;
}

// use as GQ in vcf
/*
0 proper 1 short [4 5] [8 9] (clip) [12 13] [16 17] (disc, unmap)
*/
int getSumOfMeiReads( vector<int> & count_table )
{
	vector<int> asMei = {4,5,8,9,12,13,16,17};
	int sum = 0;
	for( vector<int>::iterator it = asMei.begin(); it != asMei.end(); it++ )
		sum += count_table[*it];
	return sum;
}


void ExecuteCmd( std::string & cmd )
{
	int cmd_status = system(cmd.c_str());
	if (cmd_status != 0) {
		std::cerr << "ERROR: Fail to run: " << cmd << std::endl;
		std::cerr << "    Exit " << cmd_status << std::endl;
		exit(1);
	}
}


void GenerateDoneFile( std::string & work_dir, const char * prefix )
{
	string touch_cmd = std::string("touch ") + work_dir;
	if (work_dir[work_dir.size()-1] != '/')
		touch_cmd += "/";	
	touch_cmd += std::string(prefix) + ".Done";
	ExecuteCmd(touch_cmd);
}

bool ExistDoneFile( std::string & work_dir, const char * prefix )
{
	string fname = work_dir;
	if (work_dir[work_dir.size()-1] != '/')	
		fname += "/";
	fname += std::string(prefix) + ".Done";
	std::ifstream dfile;
	dfile.open(fname.c_str());
	if ( dfile.good() ) {
		dfile.close();
		std::cout << "Warning: exist " << fname << ", skip related steps!" << std::endl;
		return 1;
	}
	return 0;
}

// support: bwa
//			bwa-mem
string GetRemapCmd( string & full_mapper, std::string fastq_prefix, std::string ref_fasta, std::string remapSam )
{
	string mapper_name = GetFileBaseName( full_mapper );
	string path_name = GetFileDirName( full_mapper );
	string remap_cmd;
	if ( mapper_name.compare("bwa-aln") == 0 ) { // old bwa
		string mapper = path_name + "bwa";
		remap_cmd = mapper + " aln " + ref_fasta + fastq_prefix + "_1.fastq > " + remapSam + ".aln_sa1.sai; ";
		remap_cmd += mapper + " aln " + ref_fasta + fastq_prefix + "_2.fastq > " + remapSam + ".aln_sa2.sai; ";
		remap_cmd += mapper + " sampe " + ref_fasta + remapSam + ".aln_sa1.sai " + remapSam + ".aln_sa2.sai " + fastq_prefix + "_1.fastq " + fastq_prefix + "_2.fastq > " + remapSam;
	}
	else if ( mapper_name.compare( "bwa-mem" ) == 0 ) { // baw mem
		string mapper = path_name + "bwa mem";
		remap_cmd = " " + ref_fasta + fastq_prefix + "_1.fastq " + fastq_prefix + "_2.fastq > " + remapSam;
	}
	else {
		std::cerr << "ERROR: can't find command for the mapper: " << mapper_name << std::endl;
		exit(1);
	}
	return remap_cmd;
}

// get base name from full name
string GetFileBaseName( string & full_name )
{
	string base_name;
	int str_length = full_name.size();
	int slash_position = -1;
	for( int i = str_length - 2; i >= 0; i-- ) { // in case slash is the last char
		if ( full_name[i] == '/' ) {
			slash_position = i;
			break;
		}
	}
	if ( slash_position == -1 )
		base_name = full_name;
	else
		base_name = full_name.substr( slash_position + 1 );
	return base_name;
}

// get dir name ( end with slash)
string GetFileDirName( string & full_name )
{
	string dir_name;
	int str_length = full_name.size();
	int slash_position = -1;
	for( int i = str_length - 2; i >= 0; i-- ) {
		if ( full_name[i] == '/' ) {
			slash_position = i;
			break;
		}
	}
	if ( slash_position == -1 )
		dir_name = string("");
	else 
		dir_name = full_name.substr(0, slash_position + 1);
	return dir_name;	
}



















