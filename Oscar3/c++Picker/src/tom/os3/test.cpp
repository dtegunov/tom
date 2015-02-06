#include <iostream>
#include <assert.h>
#include <complex>
#include <sstream>
#include <time.h>
//-------------------------
#include "tom/os3/os3_volume.hpp"
#include "tom/os3/os3_functions.hpp"
#include "tom/os3/os3_types.hpp"
#include "tom/os3/os3_results.hpp"
#include "tom/os3/os3_options.hpp"
#include "tom/os3/os3_structures.hpp"
#include "tom/os3/os3_picklist.hpp"
#include "tom/os3/os3_picker.hpp"
//-------------------------
#include "tom/core/volume.hpp"
#include "tom/core/fftw_plan.hpp"
#include "tom/core/volume_fcn.hpp"	
#include "tom/core/io.h"
/*
std::string strI("/fs/home/hrabe/develop/images/26S/26S.em");
std::string strT("/fs/home/hrabe/develop/images/26S/26S.em");
std::string strM("/fs/home/hrabe/develop/images/26S/mask26S.em");
std::string strR("/fs/home/hrabe/develop/images/26S/res");
*/

std::string strI("/fs/home/hrabe/develop/images/26S_2/26S_106.em");
std::string strT("/fs/home/hrabe/develop/images/26S_2/26S.em");
std::string strM("/fs/home/hrabe/develop/images/26S_2/mask.em");
std::string strR("/fs/home/hrabe/develop/images/26S_2/res");

/*
std::string strI("/fs/home/hrabe/develop/images/20S/20S_3_2bin.em");
std::string strT("/fs/home/hrabe/develop/images/20S/20S_2bin.em");
std::string strM("/fs/home/hrabe/develop/images/20S/mask_2bin.em");
std::string strR("/fs/home/hrabe/develop/images/20S/res");
*/


/*
std::string strI("/fs/home/hrabe/develop/images/rib3d.em");
std::string strT("/fs/home/hrabe/develop/images/rib3d.em");
std::string strM("/fs/home/hrabe/develop/images/mask3d.em");
std::string strR("/fs/home/hrabe/develop/images/res3d.em");
*/


/*
std::string strI("/Users/thomas/workspace/images/20S_3_2bin.em");
std::string strT("/Users/thomas/workspace/images/20S_2bin.em");
std::string strM("/Users/thomas/workspace/images/mask_2bin.em");
std::string strR("/Users/thomas/workspace/images/res.em");
*/

int testSphereMask(){

	tom::Volume<TFLOAT> mask(32,32,1,NULL,NULL);

	tom::init_spheremask(mask,(float)10.0,(float).0,(float).0);
	
	mask.write_to_em("/fs/home/hrabe/develop/images/sphMask.em",NULL);

	return 0;
}

int testPSR(){

	tom::os3_volume<TFLOAT> image(strI,NULL,NULL,NULL);
	tom::os3_volume<TFLOAT> pattern(strT,NULL,NULL,NULL);
	tom::os3_volume<TFLOAT> pattern2(strT,NULL,NULL,NULL);
	tom::os3_volume<TFLOAT> mask(strM,NULL,NULL,NULL);
	tom::os3_flcf(&image,&pattern,&mask);
	tom::os3_volume<TFLOAT> flcf(&pattern);
	
	tom::os3_psr(&flcf,&pattern2);

	std::cout << pattern2.getSizeX() << std::endl;
	
	pattern2.writeToEm(std::string("/Users/thomas/workspace/images/psr.em"));
	
	return 0;
}


int testInnerSum(){

	tom::os3_volume<TFLOAT> mask(strM,NULL,NULL,NULL);
	std::cout << mask.getInnerSum() << std::endl;

	return 0;
}

int testUpdate(){

	tom::os3_results<TFLOAT> res;
	tom::Volume<TFLOAT> v(32,32,32,NULL,NULL);

	v.setValues(0);
	res.update(&v,&v,&v,0);
	
	v.setValues(1);
	res.update(&v,&v,&v,1);

	res.writeToEm(std::string("/fs/home/hrabe/develop/images/"));

	return 0;
}

int testCopy(){

	tom::os3_volume<TFLOAT> image(strT,NULL,NULL,NULL);
	tom::os3_volume<TFLOAT> mask(strM,NULL,NULL,NULL);
	
	image.calculateStatisticVolumes(&mask);

	tom::os3_volume<TFLOAT> i(&image);

	return 0;
}

int testPlanCopy(){

	tom::Volume<TFLOAT> v1(10,10,1,NULL,NULL);	
	tom::Volume<std::complex<TFLOAT> > v2(10,6,1,NULL,NULL);

	tom::fftw::Plan<TFLOAT> f1(v1,v2,FFTW_MEASURE);
	tom::fftw::Plan<TFLOAT> f2(f1);

	return 0;
}


int testAngleList(){


	tom::angleTupel start,inc,end;

	start.phi = 0;
	start.psi = 0;
	start.theta = 0;

	inc.phi = 10;
	inc.psi = 10;
	inc.theta = 10;

	end.phi = 9;
	end.psi = 9;
	end.theta = 9;
										
	std::vector<tom::angleTupel> list = tom::generateAngleList(start,inc,end);

	std::cout << list.size() << std::endl;
	tom::angleTupel t = (list.back());
	double angles[3];
	tom::angleTupelTOarray(t,&angles[0]);
	std::cout << "(" << angles[0] << "," << angles[1] << "," << angles[2] << ")" << std::endl;
	
	return 0;
}



int testSubregionWrite(){


	tom::Volume<TDOUBLE> s(20,20,20,NULL,NULL);

	s.setValues(1);
	uint32_t pos[3] = {20,20,20};

	s.write_to_em(std::string("/fs/home/hrabe/develop/images/sub.em"),NULL,&pos[0]);

	return 0;
}

int testOptions(){

	tom::os3_options opt;

	opt.setVolumeDir(std::string("/fs/home/hrabe/develop/images/20S/volumes"));
	opt.setPatternDir(std::string("/fs/home/hrabe/develop/images/20S/pattern"));
	opt.setPsfDir(std::string("/fs/home/hrabe/develop/images/20S/psf"));
	opt.setMaskDir(std::string("/fs/home/hrabe/develop/images/20S/mask"));

	opt.readDirectories();
	opt.setVolumeProperties(4,0,0,0,0,0);

	std::vector<tom::os3_job> jobList = opt.createJobList();

	return 0;

}


int testRead(){

	uint32_t s[6];
	s[0] = 0;
	s[1] = 480;
	s[2] = 0;
	s[3] = 511;
	s[4] = 512+32;
	s[5]  = 1;

	tom::os3_volume<TFLOAT> v(std::string("/fs/home/hrabe/develop/images/20S/volumes/20S_3_2bin.em"),&s[0],NULL,NULL);

	return 0;
}

int testPicklist(){

	tom::os3_job job;
	job.volumePath= std::string("/fs/home/hrabe/develop/images/20S/volumes/20S_3_2bin.em");
	job.resultPath = std::string("/fs/home/hrabe/develop/images/20S/res/20S_3_2bin.em");
	job.patternPath = std::string("/fs/home/hrabe/develop/images/20S/pattern/20S_2bin.em");
	job.jobType = __TOM_OS3_PICKLIST_JOB__;
	job.weight[0] = 1;
	job.weight[1] = 0;
	job.weight[2] = 0;
	job.numberParticles = 50;

	tom::angleTupel start,inc,end;
    start.phi = 0;
    start.psi = 0;
    start.theta = 0;
        inc.phi = 10;
        inc.psi = 10;
        inc.theta = 10;
        end.phi = 359;
        end.psi = 9;
        end.theta = 9;
										
        job.angleList = tom::generateAngleList(start,inc,end);


 		tom::os3_picklist<TFLOAT> picklist;
		
  		tom::os3_createPicklist<TFLOAT>(picklist,job,false);
		
		picklist.savePicklist(std::string("/fs/home/hrabe/develop/images/20S/list.txt"));
		std::cout << "save" << std::endl;
        picklist.loadPicklist(std::string("/fs/home/hrabe/develop/images/20S/list.txt"),false);
		std::cout << "load1" << std::endl;
		tom::os3_picklist<TFLOAT> p2;
		
		p2.loadPicklist(std::string("/fs/home/hrabe/develop/images/20S/list.txt"),false);
		std::cout << "load2" << std::endl;
		picklist.push_back(p2);
		std::cout << "merge" << std::endl;
		picklist.sortPicklist();
		std::cout << "sort" << std::endl;
		picklist.savePicklist(std::string("/fs/home/hrabe/develop/images/20S/list2.txt"));
	return 0;
}

int testStatics(){

	tom::os3_options::createDirectories(std::string("/fs/home/hrabe/develop/images/test"));

	return 0;
}

int testOptionsRead(){

	tom::os3_options options;
	options.readFromFile(std::string("/fs/home/hrabe/develop/images/test/options/options.txt"));

	std::cout << options.getVolumeDir() << std::endl;
	std::cout << options.getPatternDir() << std::endl;
	std::cout << options.getMaskDir() << std::endl;
	std::cout << options.getResultDir() << std::endl;
	std::cout << options.getPsfDir() << std::endl;
	
	
	std::vector<tom::os3_job> jobs = options.createJobList();
	
	std::cout << jobs.size() << std::endl;
	std::cout << jobs[1].angleList.size() << std::endl;
	return 0;
}

int test3dSplit(){

	tom::os3_options opt;

	tom_io_em_header volH;
	tom_io_em_read_header("/fs/home/hrabe/develop/images/3d/test/vol.em", NULL, &volH,NULL);	
	tom_io_em_header patH;
	tom_io_em_read_header("/fs/home/hrabe/develop/images/3d/pattern/rib.em", NULL, &patH,NULL);	
	
	
	tom::st_idx subVolumeSize;
	subVolumeSize.x = 64;
	subVolumeSize.y = 64;
	subVolumeSize.z = 40;

	std::pair<std::vector<tom::subVolumeCoordinates>,std::vector<tom::subVolumeCoordinates> > p;
	p = tom::split3d(&volH,&patH,subVolumeSize);

	std::vector<tom::subVolumeCoordinates> sub = p.first;
	std::vector<tom::subVolumeCoordinates> wri = p.second;
	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;
	std::cout << "-----------------------------------------------------------------" << std::endl;
	for(std::size_t i=0;i < sub.size();i++){
 		std::cout  << sub[i].startX << "," << sub[i].startY << "," << sub[i].startZ << "," << sub[i].noVoxelsX << "," << sub[i].noVoxelsY << "," << sub[i].noVoxelsZ<< std::endl;
 		std::cout << wri[i].startX << "," << wri[i].startY << "," << wri[i].startZ << "," << wri[i].noVoxelsX << "," << wri[i].noVoxelsY << "," << wri[i].noVoxelsZ << std::endl<< std::endl;
	}

	std::cout << "number subvolumes: " << wri.size() << std::endl;

	return 0;

}

int testPSF(){

	tom::os3_volume<TFLOAT> rib("/fs/home/hrabe/develop/images/3d/pattern/rib.em");
	tom::os3_volume<TFLOAT> psf("/fs/home/hrabe/develop/images/3d/psf/psf_rib.em");

	std::cout << rib.getSizeX() << rib.getSizeY() << rib.getSizeZ() << std::endl;
	std::cout << psf.getSizeX() << psf.getSizeY() << psf.getSizeZ() << std::endl;
	
	tom::os3_apply_weight_function(rib,psf);

	rib.writeToEm(std::string("/fs/home/hrabe/develop/images/3d/ribPSF.em"));

	return 0;
}

int main(int argc, char** argv){

	typedef float TFLOAT;
	//typedef double TDOUBLE;
	tom::fftw::setFcnMemFFTW<TFLOAT>();
	
        tom::os3_pick<TFLOAT> p;
        p.volumePath =";";
        p.patternPath=";";
        tom::os3_pick<TFLOAT> p1(p);
        
// 	testPSF();

 	//start3dPicker();
	//startPicker();
	//testOptions();
	//testPicklist();
	//test3dSplit();
// 	testStatics();
testOptionsRead();
	return 0;
}






























int startPicker(){
	time_t rawtime;
  	struct tm * timeinfo;
  	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	std::cout <<  "Current local time and date: " <<  asctime (timeinfo) << std::endl;
 
	tom::os3_options opt;

	opt.setVolumeDir(std::string("/fs/home/hrabe/develop/images/20S/volumes"));
	opt.setPatternDir(std::string("/fs/home/hrabe/develop/images/20S/pattern"));
	opt.setPsfDir(std::string("/fs/home/hrabe/develop/images/20S/psf"));
	opt.setMaskDir(std::string("/fs/home/hrabe/develop/images/20S/mask"));
	opt.setResultDir(std::string("/fs/home/hrabe/develop/images/20S/res/"));
	opt.readDirectories();
	opt.setVolumeProperties(2,0,0,0,0,0);

	tom::angleTupel start,inc,end;
	start.phi = 0;
	start.psi = 0;
	start.theta = 0;
	inc.phi = 10;
	inc.psi=10;
	inc.theta = 10;
	end.phi = 358;
	end.psi = 9;
	end.theta = 9;
	opt.setAngles(start,inc,end);
	std::vector<tom::os3_job> jobList = opt.createJobList();
	std::vector<tom::os3_job>::iterator iterator;
	
	std::vector<tom::os3_job> finishedJobList;
	
	std::cout << jobList.size() << std::endl;

    for( iterator = jobList.begin(); iterator != jobList.end(); iterator++){
		//(*iterator).printJob();
		tom::os3_picker<TFLOAT>(*iterator,true);
		finishedJobList.push_back(*iterator);

		
		tom::os3_job* newJob = tom::checkForFinishedVolumes(finishedJobList);
// 		std::auto_ptr<tom::os3_job> autJob(newJob);

		if(newJob->jobType == __TOM_OS3_PICKLIST_JOB__){
			//insert the stack job into the job queue
			//set the correlation filter settings of the stack job
			//setJobFilterWeight(newJob,opt.getFilterWeight());
			
			jobList.insert(iterator+1,*newJob);
		}
	}
	

	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	std::cout << std::endl <<  "Current local time and date: " <<  asctime (timeinfo) << std::endl;

	return 0;
	
}


int start3dPicker(){
	time_t rawtime;
  	struct tm * timeinfo;
  	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	std::cout <<  "Current local time and date: " <<  asctime (timeinfo) << std::endl;
 
	tom::os3_options opt;

	opt.setVolumeDir(std::string("/fs/home/hrabe/develop/images/3d/test"));
	opt.setPatternDir(std::string("/fs/home/hrabe/develop/images/3d/pattern"));
	opt.setPsfDir(std::string("/fs/home/hrabe/develop/images/3d/psf"));
	opt.setMaskDir(std::string("/fs/home/hrabe/develop/images/3d/mask"));
	opt.setResultDir(std::string("/fs/home/hrabe/develop/images/3d/res/"));
	opt.readDirectories();
	opt.setVolumeProperties(64,64,40,0,0,0);

	tom::angleTupel start,inc,end;
	start.phi = 0;
	start.psi = 0;
	start.theta = 0;
	inc.phi = 10;
	inc.psi=10;
	inc.theta = 10;
	end.phi = 9;
	end.psi = 9;
	end.theta = 9;
	opt.setAngles(start,inc,end);
	std::vector<tom::os3_job> jobList = opt.createJobList();
	std::vector<tom::os3_job>::iterator iterator;
	
	std::vector<tom::os3_job> finishedJobList;
	
	std::cout << jobList.size() << std::endl;

    for( iterator = jobList.begin(); iterator != jobList.end(); iterator++){
		//(*iterator).printJob();
		tom::os3_picker<TFLOAT>(*iterator,false);
		finishedJobList.push_back(*iterator);

		
		tom::os3_job* newJob = tom::checkForFinishedVolumes(finishedJobList);
// 		std::auto_ptr<tom::os3_job> autJob(newJob);

		if(newJob->jobType == __TOM_OS3_PICKLIST_JOB__){
			//insert the stack job into the job queue
			//set the correlation filter settings of the stack job
			//setJobFilterWeight(newJob,opt.getFilterWeight());
			
// 			jobList.insert(iterator+1,*newJob);
		}
	}
	
	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	std::cout << std::endl <<  "Current local time and date: " <<  asctime (timeinfo) << std::endl;

	return 0;
	
}










