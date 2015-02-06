/****************************************************************************//**
 * \file tom_os3_structures.hpp
 * \brief The header file defining important structures
 * \author  Thomas Hrabe
 * \version 0.1
 * \date    08.01.2008
 *******************************************************************************/ 
#ifndef __INCLUDE_TOM_OS3_STRUCTURES_HPP__
#define __INCLUDE_TOM_OS3_STRUCTURES_HPP__

#define PI 3.14159265358979323846264338327950288419716939937510

#include "tom/os3/os3_types.hpp"

#include "tom/core/volume_fcn.hpp"

#include <string>
#include <vector>

namespace tom{
	
	

	typedef struct angleTupel{
		double phi;
		double psi;
		double theta;
	
		angleTupel(){};
		angleTupel(const tom::angleTupel& original){
			this->phi =original.phi;
			this->psi =original.psi;
			this->theta =original.theta;
		};

	} angleTupel;
	void angleTupelTOarray(const tom::angleTupel tupel,double* array);
	
	template<typename T>
	struct os3_pick {
		std::string volumePath;
		std::string patternPath;
		st_idx position;
		tom::angleTupel angles;
		
		T value;
		
		os3_pick(){}
		os3_pick(std::string volume,st_idx pos,tom::angleTupel angles,T value){
			this->volumePath = volume;
			this->position = position;
			this->angles = angles;
			this->value = value;
		}
		os3_pick(const tom::os3_pick<T>& original){
			this->volumePath = original.volumePath;
                        this->patternPath = original.patternPath;
			this->position.x = original.position.x;
                        this->position.y = original.position.y;
                        this->position.z = original.position.z;
			this->angles = original.angles;
			this->value = original.value;
		}

		bool operator<(const os3_pick& other) const{	
			//sort in ascending order!
			return this->value > other.value;
		}

	};

	typedef struct subVolumeCoordinates{

		uint32_t startX;
		uint32_t startY;
		uint32_t startZ;

		uint32_t noVoxelsX;
		uint32_t noVoxelsY;
		uint32_t noVoxelsZ;

	}subVolumeCoordinates;

	typedef struct os3_job{
			int jobType;

			std::string volumePath;
			std::string patternPath;
			std::string maskPath;
			std::string psfPath;
			std::string resultPath;
			
			/*settings for the correlation job*/
			uint32_t subregion[6];
			uint32_t writeBackSubregion[6];
			uint32_t binning[3];
			uint32_t sampling[3];
			std::size_t id;
			std::size_t numberJobs;
			std::vector<tom::angleTupel> angleList;

			/*settings for the stack job*/
			double weight[3];		
			std::size_t numberParticles;

			bool filterSave;
			bool listSave;
			bool stackSave;
			/*settings for the classification job*/
			
			/*constructors*/
			os3_job(){};
			os3_job(const tom::os3_job& originalJob){
				this->jobType = originalJob.jobType;
				this->volumePath = originalJob.volumePath;
				this->patternPath = originalJob.patternPath;
				this->maskPath = originalJob.maskPath;
				this->psfPath = originalJob.psfPath;
				this->resultPath = originalJob.resultPath;
				this->id = originalJob.id;
				this->numberJobs = originalJob.numberJobs;
				this->angleList = originalJob.angleList;

				this->subregion[0] = originalJob.subregion[0];
				this->subregion[1] = originalJob.subregion[1];
				this->subregion[2] = originalJob.subregion[2];
				this->subregion[3] = originalJob.subregion[3];
				this->subregion[4] = originalJob.subregion[4];
				this->subregion[5] = originalJob.subregion[5];

				this->writeBackSubregion[0] = originalJob.writeBackSubregion[0];
				this->writeBackSubregion[1] = originalJob.writeBackSubregion[1];
				this->writeBackSubregion[2] = originalJob.writeBackSubregion[2];
				this->writeBackSubregion[3] = originalJob.writeBackSubregion[3];
				this->writeBackSubregion[4] = originalJob.writeBackSubregion[4];
				this->writeBackSubregion[5] = originalJob.writeBackSubregion[5];

				this->binning[0] = originalJob.binning[0];
				this->binning[1] = originalJob.binning[1];
				this->binning[2] = originalJob.binning[2];

				this->sampling[0] = originalJob.sampling[0];
				this->sampling[1] = originalJob.sampling[1];
				this->sampling[2] = originalJob.sampling[2];

				this->weight[0] = originalJob.weight[0];
				this->weight[1] = originalJob.weight[1];
				this->weight[2] = originalJob.weight[2];

				this->numberParticles = originalJob.numberParticles;

				this->filterSave = originalJob.filterSave;
				this->listSave = originalJob.listSave;
				this->stackSave = originalJob.stackSave;
			}

	}os3_job;


void setJobSubregion(tom::os3_job& job,tom::subVolumeCoordinates& coordinates,tom::subVolumeCoordinates& writeBackCoordinates);
void setJobSampling(tom::os3_job& job,uint32_t* sampling);
void setJobBinning(tom::os3_job& job,uint32_t* binning);
void setJobFilterWeight(tom::os3_job& job,double* weight);
void jobInfo(tom::os3_job& job);

}

template struct tom::os3_pick<TFLOAT>;
template struct tom::os3_pick<TDOUBLE>;


	inline void tom::angleTupelTOarray(const tom::angleTupel tupel,double* array){
	
		array[0] = tupel.phi;
		array[1] = tupel.psi;
		array[2] = tupel.theta;
	}

#endif
