/*
    Copyright (c) 2016, Konstantinos Kamnitsas
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the copyright holder nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/


#include "densecrf.h"
#include <cstdio>
#include <cmath>

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMetaDataObject.h"
#include "itkNiftiImageIO.h"

#include <vector>
#include <string>
#include <sstream> //for ostringstream.
#include <fstream>

typedef std::vector<std::string> vectorOfStringsType;	
typedef std::vector<vectorOfStringsType> vectorOfVectorOfStringsType;

typedef float InputPixelType;
const unsigned int  DimensionOfNii = 3;
typedef itk::Image< InputPixelType, DimensionOfNii > InputImageType;

typedef itk::NiftiImageIO ImageIOType;

    
/*
 * --------------------Loader and saver of niftis, as well the function that takes a resulting (from inference matrix/vector) and creates a nii out of it.-----------------------
 */
InputImageType::Pointer loadNiiImageItk( std::string niiFullFilepathname, ImageIOType::Pointer & niftiImageIO ) {
    //http://www.itk.org/Doxygen/html/classitk_1_1Image.html
    //http://www.itk.org/Doxygen/html/classitk_1_1ImageFileReader.html
    
    typedef itk::ImageFileReader<InputImageType> ReaderType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(niiFullFilepathname);
    

    //ImageIOType::Pointer niftiImageIO = ImageIOType::New();
    reader->SetImageIO( niftiImageIO );
    
    try {
        reader->Update();
    } catch (itk::ExceptionObject& e) {
        std::cerr << e.GetDescription() << std::endl;
        exit(1);  // You can choose to do something else, of course.
    }

    //InputImageType::Pointer inputImage = reader->GetOutput();
    InputImageType::Pointer inputImage = reader->GetOutput();
    
    return inputImage;

}


int saveNiiImageItk( std::string niiFullFilepathname, InputImageType::Pointer & outputImageItkType , ImageIOType::Pointer & niftiImageIO) {
  std::cout << "Saving image to: " << niiFullFilepathname << "\n";

  typedef itk::ImageFileWriter< InputImageType >  Writer1Type;
  Writer1Type::Pointer writer1 = Writer1Type::New();

  writer1->SetInput( outputImageItkType );
  writer1->SetFileName( niiFullFilepathname );
  writer1->SetImageIO( niftiImageIO ); //seems like this is useless.

  // Execution of the writer is triggered by invoking the \code{Update()} method.
  try
    {
    writer1->Update();
    }
  catch (itk::ExceptionObject & e)
    {
    std::cerr << "exception in file writer " << std::endl;
    std::cerr << e.GetDescription() << std::endl;
    std::cerr << e.GetLocation() << std::endl;
    return 1;
    }

  return 0;
}

 
void makeNiiProbMapFromMatrixXfOfProbsAndSave(std::string folderForThePatientToSaveTheResultsIn,
											std::string stringWithThePrefixOfTheFilenameToSaveResultingProbMaps, 
											MatrixXf & probMaps1, //first row is the background class.
											int numberOfForegroundClasses, //excluding background.
											int * sizeOfImages, 
											InputImageType::Pointer imgOneOfModsItkTypeForReferenceOfResultingImage,
											ImageIOType::Pointer niftiImageIOForReferenceOfResultingImage) {
												
	int	numberOfVoxelsInImage = sizeOfImages[0]*sizeOfImages[1]*sizeOfImages[2];

    typedef InputImageType::IndexType InputImageIndexType; //should be a 3 dimensional array, one for each dimension of the image.
    typedef InputImageType::PixelType InputImagePixelType; //this will return Double, or whatever I defined at the very top when I define the InputImageType?!
    InputImageIndexType index_img;
	InputImagePixelType voxel_val;				


	for (int class_i=0; class_i <= numberOfForegroundClasses ; class_i++) { //class_i == 0 is the background, which is not counted in the numberOfForegroundClasses.
		//lets create an image with this class's probabilities. Use as a starting point/reference an image-itk-pointer passed from the modalities.
		//if it is the background, skip this, dont create an image:
		//if (class_i==0)
			//continue;
		//BEAR IN MIND THAT THIS CHANGES THE VERY ORIGINAL IMAGE IN MEMORY, SO DON'T USE IT AGAIN AFTER THIS POINT!

		for (int z=0; z<sizeOfImages[2]; z++) {
			for (int c=0; c<sizeOfImages[1]; c++) {
				for (int r=0; r<sizeOfImages[0]; r++) {
					index_img[0] = r; index_img[1] = c; index_img[2] = z;
					voxel_val = probMaps1(class_i, r + c*sizeOfImages[0] + z*sizeOfImages[0]*sizeOfImages[1]);

					imgOneOfModsItkTypeForReferenceOfResultingImage->SetPixel(index_img, voxel_val);
				}
			}
		}
		
		//now save it.
		std::ostringstream convertNumToStr;
		convertNumToStr << class_i;
		std::string niiFullFilepathname = folderForThePatientToSaveTheResultsIn + "/" + stringWithThePrefixOfTheFilenameToSaveResultingProbMaps + convertNumToStr.str() + ".nii.gz";
		saveNiiImageItk( niiFullFilepathname, imgOneOfModsItkTypeForReferenceOfResultingImage , niftiImageIOForReferenceOfResultingImage);
	}
}

void makeNiiSegmMapFromVectorXsOfLabelsAndSave(std::string folderForThePatientToSaveTheResultsIn,
											std::string stringWithThePrefixOfTheFilenameToSaveResultingSegmMap,
											VectorXs & map1, //first row is the background class.
											int * sizeOfImages, 
											InputImageType::Pointer imgOneOfModsItkTypeForReferenceOfResultingImage,
											ImageIOType::Pointer niftiImageIOForReferenceOfResultingImage) {
												
	int	numberOfVoxelsInImage = sizeOfImages[0]*sizeOfImages[1]*sizeOfImages[2];

    typedef InputImageType::IndexType InputImageIndexType; //should be a 3 dimensional array, one for each dimension of the image.
    typedef InputImageType::PixelType InputImagePixelType; //this will return Double, or whatever I defined at the very top when I define the InputImageType?!
    InputImageIndexType index_img;
	InputImagePixelType voxel_val;				
			  		
	for (int z=0; z<sizeOfImages[2]; z++) {
			for (int c=0; c<sizeOfImages[1]; c++) {
				for (int r=0; r<sizeOfImages[0]; r++) {
					index_img[0] = r; index_img[1] = c; index_img[2] = z;
					voxel_val = map1( r + c*sizeOfImages[0] + z*sizeOfImages[0]*sizeOfImages[1] );

					imgOneOfModsItkTypeForReferenceOfResultingImage->SetPixel(index_img, voxel_val);
				}
			}
	}
		
	//now save it.
	std::string niiFullFilepathname = folderForThePatientToSaveTheResultsIn + "/" + stringWithThePrefixOfTheFilenameToSaveResultingSegmMap + ".nii.gz";
	saveNiiImageItk( niiFullFilepathname, imgOneOfModsItkTypeForReferenceOfResultingImage , niftiImageIOForReferenceOfResultingImage);
	
}
											
											
//Helper function. The algorithm works in this 0-250 range in order to be efficient.			
float transFTo0255(float f, float minIntensity, float maxIntensity) {
    double trans = ( 255.0 / (maxIntensity - (minIntensity) ) ) * ( f - ( minIntensity ));
    if (trans < 0)
        trans = 0;
    else if (trans > 254)
        trans = 254;
    return trans;
}

// Loads a nii and make the matrix that the rest of the code needs.
unsigned char * readNiiModalitiesOfPatientToNeededFormat( vectorOfStringsType patientModalitiesFullPathNames,
															int numModalities,
															int * sizeOfImages,
															float minIntensity,
															float maxIntensity,
															InputImageType::Pointer & imgOneOfModsItkTypeForReferenceOfResultingImage,
															ImageIOType::Pointer & niftiImageIOForReferenceOfResultingImage) {
    std::vector<InputImageType::Pointer> vectorWithAllModalityImagesItkTypeForThisPatient;
    
	for (int modality_i = 0 ; modality_i < numModalities ; modality_i ++) {
		std::cout<<"Loading modality: " <<patientModalitiesFullPathNames[modality_i]<<"\n";
		ImageIOType::Pointer niftiImageIO = ImageIOType::New();
		InputImageType::Pointer inputImageFromItk = loadNiiImageItk( patientModalitiesFullPathNames[modality_i], niftiImageIO);
		
		vectorWithAllModalityImagesItkTypeForThisPatient.push_back(inputImageFromItk);
	    //Lets get the size of the image and check that is correct with what is given as configuration:
		typedef InputImageType::SizeType InputImageSize;
		InputImageSize img_size = inputImageFromItk->GetLargestPossibleRegion().GetSize(); //Returns the size of the whole 3D image.
		int R = img_size[0]; //eg 230, 230, 154
		int C = img_size[1];
		int Z = img_size[2];
		if ( R != sizeOfImages[0] || C != sizeOfImages[1] || Z != sizeOfImages[2] ) {
			std::cout<<"ERROR! Size of the image read is not the same as given as a configuration! Exiting!\n";
			std::cout<<"ERROR for image:"<<patientModalitiesFullPathNames[modality_i]<<"\n";
			exit(1);
		}
		
		
		//Pass back 1 image and 1 imageIOType, to use as reference when saving the resulting probability-maps.
		if (modality_i==0) {
			imgOneOfModsItkTypeForReferenceOfResultingImage = inputImageFromItk;
			niftiImageIOForReferenceOfResultingImage = niftiImageIO;
		}
	}

	std::cout<<"All the modalities were successfully read for patient.\n";

    
    //Create the array that will hold the image-modalities data, as wanted by the algorithm. The channels must be in the same array.
    
    int numberOfVoxelsInImageMultipliedByModalities = sizeOfImages[0]*sizeOfImages[1]*sizeOfImages[2]*numModalities;
    unsigned char * arrayWithAllModalitiesOfPatient = new unsigned char[numberOfVoxelsInImageMultipliedByModalities];
    
    //Defining the type of index and the type of value that the pixels have.
    typedef InputImageType::IndexType InputImageIndexType; //should be a 3 dimensional array, one for each dimension of the image.
    typedef InputImageType::PixelType InputImagePixelType; //this will return Double, or whatever I defined at the very top when I define the InputImageType?!
    
    //define the variables
    InputImageIndexType index_img;
	InputImagePixelType targetPixel_val;
	
	for (int pixelIndexInResultingArray_i = 0 ; pixelIndexInResultingArray_i < numberOfVoxelsInImageMultipliedByModalities ; pixelIndexInResultingArray_i++) {
        //Calculate the coordinates of the voxel that I want to touch, in the actual images.
		int actualR=-1, actualC=-1, actualZ=-1, actualModality_i=-1; //these are all indeces. Eg actualModality_i = 1 corresponds to the second available modality.
		
		//they change from faster to slower: channel, R, C, Z.
		actualModality_i = pixelIndexInResultingArray_i % numModalities ;
		actualR = (pixelIndexInResultingArray_i / numModalities) % sizeOfImages[0] ;
		actualC = (pixelIndexInResultingArray_i / ( numModalities * sizeOfImages[0] ) ) % sizeOfImages[1] ;
		actualZ = (pixelIndexInResultingArray_i / ( numModalities * sizeOfImages[0] * sizeOfImages[1] ) ) % sizeOfImages[2] ;
		
		index_img[0] = actualR;
		index_img[1] = actualC;
		index_img[2] = actualZ;
		
		targetPixel_val = vectorWithAllModalityImagesItkTypeForThisPatient[actualModality_i]->GetPixel(index_img);

        //Here I need to normalize the pixel values to the 0-255 regime...

        //float -3 -> +3 normalize to 0 -> 255
        unsigned char targetPixel_val_0to255 = (unsigned char) transFTo0255(targetPixel_val, minIntensity, maxIntensity);

        arrayWithAllModalitiesOfPatient[pixelIndexInResultingArray_i] = targetPixel_val_0to255;

	}


    std::cout<<"Done with loading the nifti modalities.\n";

    return arrayWithAllModalitiesOfPatient;
}



	


//Loads the nii with the unary potentials and makes a matrix that the rest of the original code needs.
MatrixXf readNiiProbMapsOfPatientAndMakeUnaryEnergyMatrix( vectorOfStringsType patientProbabilityMapsFullPathNames, int numberOfForegroundClasses, int * sizeOfImages) {
    std::vector<InputImageType::Pointer> vectorWithAllProbabilityImagesItkTypeForThisPatient;

    for (int class_i = 0 ; class_i < numberOfForegroundClasses ; class_i ++) {
		ImageIOType::Pointer niftiImageIO = ImageIOType::New();
        InputImageType::Pointer inputImageFromItk = loadNiiImageItk( patientProbabilityMapsFullPathNames[class_i], niftiImageIO );
        vectorWithAllProbabilityImagesItkTypeForThisPatient.push_back(inputImageFromItk);
        //Lets get the size of the image and check that is correct with what is given as configuration:
        typedef InputImageType::SizeType InputImageSize;
        InputImageSize img_size = inputImageFromItk->GetLargestPossibleRegion().GetSize(); //Returns the size of the whole 3D image.
        int R = img_size[0]; //eg 230, 230, 154
        int C = img_size[1];
        int Z = img_size[2];
        if ( R != sizeOfImages[0] || C != sizeOfImages[1] || Z != sizeOfImages[2] ) {
            std::cout<<"ERROR! Size of the image read is not the same as given as a configuration! Exiting!\n";
            std::cout<<"ERROR for image:"<<patientProbabilityMapsFullPathNames[class_i]<<"\n";
            exit(1);
        }
    }
    std::cout<<"Alright, all the modalities were successfully read for patient.\n";


    //Create the array that will hold the image-probability data, as wanted by the algorithm. Columns are the classes (background included), rows are the pixels (rcz, as with images).

    int numberOfVoxelsInImage = sizeOfImages[0]*sizeOfImages[1]*sizeOfImages[2];

    MatrixXf matrixWithUnaries( numberOfForegroundClasses+1, numberOfVoxelsInImage );


    //Defining the type of index and the type of value that the pixels have.
    typedef InputImageType::IndexType InputImageIndexType; //should be a 3 dimensional array, one for each dimension of the image.
    typedef InputImageType::PixelType InputImagePixelType; //this will return Double, or whatever I defined at the very top when I define the InputImageType?!

    //define the variables
    InputImageIndexType index_img;
    InputImagePixelType targetPixel_val;


    for (int pixelIndexInResultingArray_i = 0 ; pixelIndexInResultingArray_i < numberOfVoxelsInImage ; pixelIndexInResultingArray_i++) {
        //Calculate the coordinates of the voxel that I want to touch, in the actual images.
        int actualR=-1, actualC=-1, actualZ=-1; //these are all indeces. Eg actualModality_i = 1 corresponds to the second available modality.

        //they change from faster to slower: channel, R, C, Z.
        actualR = pixelIndexInResultingArray_i % sizeOfImages[0] ;
        actualC = (pixelIndexInResultingArray_i / sizeOfImages[0]) % sizeOfImages[1] ;
        actualZ = (pixelIndexInResultingArray_i / ( sizeOfImages[0] * sizeOfImages[1] ) ) % sizeOfImages[2] ;

        index_img[0] = actualR;
        index_img[1] = actualC;
        index_img[2] = actualZ;

        //std::cout<<index_img[0]<<" "<<index_img[1]<<" "<<index_img[2]<<" "<<actualModality_i<<"\n";
        float totalProbabilityOfClassesWithoutBackgr = 0;
        for (int class_i=0 ; class_i < numberOfForegroundClasses ; class_i++) { //numberOfForegroundClasses is without counting the background.
            targetPixel_val = vectorWithAllProbabilityImagesItkTypeForThisPatient[class_i]->GetPixel(index_img);
            matrixWithUnaries(class_i+1, pixelIndexInResultingArray_i) = -log(targetPixel_val); //class_i+1 because I am keeping the first slot for the background class! :)
            totalProbabilityOfClassesWithoutBackgr += targetPixel_val;
        }
        //Now lets add the rest of the probability into the background class, which will be in the very last row of the matrix...
        float backgrProb = 1.0 - totalProbabilityOfClassesWithoutBackgr;
        if (backgrProb<0)
            std::cout<<"WARN: BACKGROUND PROBABILITY WAS CALCULATED NEGATIVE! ZEROING IT OUT, BUT CHECK IT!!!!";
        backgrProb = backgrProb<0? 0 : backgrProb;
        matrixWithUnaries(0, pixelIndexInResultingArray_i) = -log(backgrProb);

    }

    std::cout<<"End of the function that returns the unaries map.\n";
    return matrixWithUnaries;

}


/*
 * -------------------- Class that holds the parameters of the CRF-----------------------
 */

class ParametersOfRun {
	
	public :
	int MaxIterations;
	float PosRStd;
	float PosCStd;
	float PosZStd;
	float PosW;
	float BilateralRStd; 
	float BilateralCStd; 
	float BilateralZStd; 
	float * BilateralModsStds;
	int numberOfModalities;
	float BilateralW; 
	
	ParametersOfRun() {
		MaxIterations = 5;
		PosRStd = 3;
		PosCStd = 3;
		PosZStd = 3;
		PosW = 3;
		BilateralRStd = 5;
		BilateralCStd = 5;
		BilateralZStd = 5;
		BilateralW = 5;
		numberOfModalities = 0;
	}
	~ParametersOfRun() {
		delete [] BilateralModsStds;
	}
	
	void setNumberOfModalitiesAndInitialiseBilateralOfModalities(int numbOfModalities) {
		numberOfModalities = numbOfModalities;
		BilateralModsStds = new float[numberOfModalities];
		for (int mod_i=0; mod_i < numberOfModalities; mod_i++) {
			BilateralModsStds[mod_i] = 3;
		}
	}
	
	void printParameters() {
		std::cout << "========== Parameters for the CRF: ============" << std::endl;
		std::cout << "MaxIterations: " << MaxIterations << std::endl;
		
		std::cout << "PosRStd: "  << PosRStd << std::endl;
		std::cout << "PosCStd: "  << PosCStd << std::endl;
		std::cout << "PosZStd: "  << PosZStd << std::endl;
		std::cout << "PosW: "     << PosW    << std::endl;
		
		std::cout << "BilateralRStd: " << BilateralRStd << std::endl;
		std::cout << "BilateralCStd: " << BilateralCStd << std::endl;
		std::cout << "BilateralZStd: " << BilateralZStd << std::endl;
		std::cout << "BilateralW: "     << BilateralW    << std::endl;
		
		std::cout << "numberOfModalities: " << numberOfModalities << std::endl;
		for (int mod_i=0 ; mod_i<numberOfModalities ; mod_i++) {
			std::cout << "BilateralModsStds[" << mod_i <<"]: " << BilateralModsStds[mod_i] << std::endl;
		}

	}
};


void addPairwiseBilateralToCrfForCorrespondingNumberOfModalities( DenseCRF3D & crf3d, 
									ParametersOfRun & parametersOfRun,
									unsigned char * arrayWithMyImWithAllModalities) {
	if (parametersOfRun.numberOfModalities==1) {
		std::cout<<"Using addPairwiseBilateral for 1 Modality\n";
		crf3d.addPairwiseBilateral1Mod( parametersOfRun.BilateralRStd,
										parametersOfRun.BilateralCStd,
										parametersOfRun.BilateralZStd,
										parametersOfRun.BilateralModsStds[0],
										arrayWithMyImWithAllModalities,
										new PottsCompatibility( parametersOfRun.BilateralW ) );
	} else if (parametersOfRun.numberOfModalities==2) {
		std::cout<<"Using addPairwiseBilateral for 2 Modalities\n";
		crf3d.addPairwiseBilateral2Mod( parametersOfRun.BilateralRStd,
										parametersOfRun.BilateralCStd,
										parametersOfRun.BilateralZStd,
										parametersOfRun.BilateralModsStds[0],
										parametersOfRun.BilateralModsStds[1],
										arrayWithMyImWithAllModalities,
										new PottsCompatibility( parametersOfRun.BilateralW ) );
	} else if (parametersOfRun.numberOfModalities==3) {
		std::cout<<"Using addPairwiseBilateral for 3 Modalities\n";
		crf3d.addPairwiseBilateral3Mod( parametersOfRun.BilateralRStd,
										parametersOfRun.BilateralCStd,
										parametersOfRun.BilateralZStd,
										parametersOfRun.BilateralModsStds[0],
										parametersOfRun.BilateralModsStds[1],
										parametersOfRun.BilateralModsStds[2],
										arrayWithMyImWithAllModalities,
										new PottsCompatibility( parametersOfRun.BilateralW ) );	
	} else if (parametersOfRun.numberOfModalities==4) {
		std::cout<<"Using addPairwiseBilateral for 4 Modalities\n";
		crf3d.addPairwiseBilateral4Mod( parametersOfRun.BilateralRStd,
										parametersOfRun.BilateralCStd,
										parametersOfRun.BilateralZStd,
										parametersOfRun.BilateralModsStds[0],
										parametersOfRun.BilateralModsStds[1],
										parametersOfRun.BilateralModsStds[2],
										parametersOfRun.BilateralModsStds[3],
										arrayWithMyImWithAllModalities,
										new PottsCompatibility( parametersOfRun.BilateralW ) );
	} else if (parametersOfRun.numberOfModalities==5) {
		std::cout<<"Using addPairwiseBilateral for 5 Modalities\n";
		crf3d.addPairwiseBilateral5Mod( parametersOfRun.BilateralRStd,
										parametersOfRun.BilateralCStd,
										parametersOfRun.BilateralZStd,
										parametersOfRun.BilateralModsStds[0],
										parametersOfRun.BilateralModsStds[1],
										parametersOfRun.BilateralModsStds[2],
										parametersOfRun.BilateralModsStds[3],
										parametersOfRun.BilateralModsStds[4],
										arrayWithMyImWithAllModalities,
										new PottsCompatibility( parametersOfRun.BilateralW ) );
	}
									
}

/*
 * ------ Custom parser of the arguments.. I should definitely use a library for this... This is ugly and not robust. ---------
 * Takes as input a vector of strings, which are all the -options and corresponding arguments. The vector is generated by readParamsFromConfigFile() or readArgvIntoVectorOfStrings()
 */
void parseArgs(
			vectorOfStringsType & vectorWithArgs,
			vectorOfStringsType & vectorWithAllModalitiesPaths, //if you dont pass this by reference, it creates a copy of the vector and the contents that are added here are not propagated to the back vector!
			vectorOfStringsType & vectorWithAllProbabilityMapsPaths,
			std::string & stringWithTheFolderForThePatientToSaveTheResultsIn,
			std::string & stringWithThePrefixOfTheFilenameToSaveResultingSegmMap,
			std::string & stringWithThePrefixOfTheFilenameToSaveResultingProbMaps,
			
			int & numberOfModalities,
			int & numberOfForegroundClasses,
			int * sizeOfImages,
			float & minIntensity,
			float & maxIntensity,
			ParametersOfRun & parametersOfRun) {

	std::cout <<"==============================================================================\n";
	std::cout <<"=============== Parsing configuration parameters for the process =============\n";
	std::cout <<"==============================================================================\n";
				
	std::string patientMainFolder = "placeholder";
	int numberOfDimensions = 0;

	for (int arg_i=0; arg_i < vectorWithArgs.size() ; arg_i++) {
		
		if (vectorWithArgs[arg_i].compare("\n") == 0 || vectorWithArgs[arg_i].compare("") == 0 || vectorWithArgs[arg_i][0]== '#') {
			//empty line or comment, continue with next line.
			continue;
		}
		
		std::cout <<"Parsing input arguments for option:\'"<<vectorWithArgs[arg_i] <<"\'\n";
		
		if (vectorWithArgs[arg_i].compare("-patientMainFolder") == 0) {
			if 	(patientMainFolder == "shouldHaveBeenReadAlready") {
				std::cout <<"ERROR: Problem when parsing the dense_inference input arguments! The option \'-patientFolder\' Should be placed before the modalities and probability map file names! Exiting!\n";
				exit(1);
			}
			
			arg_i++;
			patientMainFolder = vectorWithArgs[arg_i];
			std::cout << "Parsed Parameter: as a main-folder for the patient: " << patientMainFolder << "\n";
		}
		else if (vectorWithArgs[arg_i].compare("-numberOfModalitiesAndFiles") == 0) {
			arg_i++;
			numberOfModalities = atoi(vectorWithArgs[arg_i].c_str());
			std::cout << "Parsed Parameter: the number of modalities to use: " << numberOfModalities << "\n";
			
			for (int mod_i=0 ; mod_i < numberOfModalities; mod_i++) {
				arg_i++;
				if (patientMainFolder == "placeholder" || patientMainFolder == "shouldHaveBeenReadAlready") {
					//Option -patientFolder with the main folder of a patient was not provided. In this case, expect to read full-filenamePaths of each modality and probability map.
					vectorWithAllModalitiesPaths.push_back(vectorWithArgs[arg_i]);
					patientMainFolder = "shouldHaveBeenReadAlready";
				}
				else {
					//Option -patientFolder with the main folder of a patient was previously provided. In this case, expect to read only filenames of modalities and prob-maps, and look for them in the main patient-folder.
					vectorWithAllModalitiesPaths.push_back(patientMainFolder + "/" + vectorWithArgs[arg_i]);
				}
				std::cout << "Parsed Parameter: the full filename+path for the input Modality-" << mod_i << " was parsed to be: " << vectorWithAllModalitiesPaths.back() << "\n";

			}
			parametersOfRun.setNumberOfModalitiesAndInitialiseBilateralOfModalities(numberOfModalities);
		}
		else if (vectorWithArgs[arg_i].compare("-numberOfForegroundClassesAndProbMapFiles") == 0) {
			arg_i++;
			numberOfForegroundClasses = atoi(vectorWithArgs[arg_i].c_str());
			std::cout << "Parsed Parameter: The number of foreground classes was parsed to be: " << numberOfForegroundClasses << "\n";
			for (int class_i=0 ; class_i < numberOfForegroundClasses; class_i++) {
				arg_i++;
				if (patientMainFolder == "placeholder" || patientMainFolder == "shouldHaveBeenReadAlready") {
					//Option -patientFolder with the main folder of a patient was not provided. In this case, expect to read full-filenamePaths of each modality and probability map.
					vectorWithAllProbabilityMapsPaths.push_back(vectorWithArgs[arg_i]);
					patientMainFolder = "shouldHaveBeenReadAlready";
				}
				else {
					//Option -patientFolder with the main folder of a patient was previously provided. In this case, expect to read only filenames of modalities and prob-maps, and look for them in the main patient-folder.
					vectorWithAllProbabilityMapsPaths.push_back(patientMainFolder + "/" + vectorWithArgs[arg_i]);
				}
				std::cout << "Parsed Parameter: the full filename+path for the input Probability-Map for Class-" << class_i << " was parsed to be: " << vectorWithAllProbabilityMapsPaths.back() << "\n";
			}
		}
		else if (vectorWithArgs[arg_i].compare("-imageDimensions") == 0) {
			//First should be given 3 if 3D or 2 if 2D image, followed by the size of the image in the X-Y(-Z) direction. 
			arg_i++;
			numberOfDimensions = atoi(vectorWithArgs[arg_i].c_str());
			std::cout << "Parsed Parameter: The image was given to be of : " << numberOfDimensions <<" Dimensions. The size of each dimension (eg R-C-Z) was parsed: ";
			for (int dim_i=0 ; dim_i < numberOfDimensions; dim_i++) {
				arg_i++;
				sizeOfImages[dim_i] = atoi(vectorWithArgs[arg_i].c_str());
				std::cout << sizeOfImages[dim_i] << " ";
			}
			std::cout << "\n";
		}
		else if (vectorWithArgs[arg_i].compare("-minMaxIntensities") == 0) {
			//First should be given 3 if 3D or 2 if 2D image, followed by the size of the image in the X-Y(-Z) direction. 
			minIntensity = atof(vectorWithArgs[arg_i+1].c_str());
			maxIntensity = atof(vectorWithArgs[arg_i+2].c_str());
			arg_i += 2;
			std::cout << "Parsed Parameter: The min and max intensities that will be processed from the channels (modalities should be scaled at the same range) were parsed to be: [" << minIntensity << "," << maxIntensity << "]\n";
		}
		else if (vectorWithArgs[arg_i].compare("-outputFolder") == 0) {
			//Output folder, where the results of the inference will be placed. The folder should already exist.
			arg_i++;
			stringWithTheFolderForThePatientToSaveTheResultsIn = vectorWithArgs[arg_i] + "/";
			std::cout << "Parsed Parameter: The output folder were results will be saved was parsed to be: " << stringWithTheFolderForThePatientToSaveTheResultsIn << "\n";
		}
		else if (vectorWithArgs[arg_i].compare("-prefixForOutputSegmentationMap") == 0) {
			arg_i++;
			stringWithThePrefixOfTheFilenameToSaveResultingSegmMap = vectorWithArgs[arg_i];
			std::cout << "Parsed Parameter: The segmentation map will be saved with the filename: " << stringWithThePrefixOfTheFilenameToSaveResultingSegmMap << "\n";
		}
		else if (vectorWithArgs[arg_i].compare("-prefixForOutputProbabilityMaps") == 0) {
			arg_i++;
			stringWithThePrefixOfTheFilenameToSaveResultingProbMaps = vectorWithArgs[arg_i];
			std::cout << "Parsed Parameter: The resulting probability maps for each class will be saved with the prefix: " << stringWithThePrefixOfTheFilenameToSaveResultingProbMaps << "\n";
		}
		//NOW THE PARAMETERS OF THE CRF:
		else if (vectorWithArgs[arg_i].compare("-pRCZandW") == 0) {
			// PosRStd, PosCStd, PosZStd, PosW
			parametersOfRun.PosRStd = atof(vectorWithArgs[arg_i+1].c_str());
			parametersOfRun.PosCStd = atof(vectorWithArgs[arg_i+2].c_str());
			parametersOfRun.PosZStd = atof(vectorWithArgs[arg_i+3].c_str());
			parametersOfRun.PosW = atof(vectorWithArgs[arg_i+4].c_str());
			arg_i += 4;
			std::cout << "Parsed Parameter: Positional-only Stds for R-C-Z and Positional-only W parsed: " << parametersOfRun.PosRStd <<
				" "<< parametersOfRun.PosCStd << " " << parametersOfRun.PosZStd << " " << parametersOfRun.PosW << "\n";
		}
		else if (vectorWithArgs[arg_i].compare("-bRCZandW") == 0) {
			// BilateralRStd, BilateralCStd, BilateralZStd, ( BilateralMod1Std, BilateralMod2Std, ... ), BilateralW
			parametersOfRun.BilateralRStd = atof(vectorWithArgs[arg_i+1].c_str());
			parametersOfRun.BilateralCStd = atof(vectorWithArgs[arg_i+2].c_str());
			parametersOfRun.BilateralZStd = atof(vectorWithArgs[arg_i+3].c_str());
			parametersOfRun.BilateralW = atof(vectorWithArgs[arg_i+4].c_str());
			arg_i+=4;
			std::cout << "Parsed Parameter: Bilateral positional Stds for R-C-Z and Bilateral W parsed: " << parametersOfRun.BilateralRStd <<
				" "<< parametersOfRun.BilateralCStd << " " << parametersOfRun.BilateralZStd << " " << parametersOfRun.BilateralW << "\n";
		}
		else if (vectorWithArgs[arg_i].compare("-bMods") == 0) {
			std::cout << "Parsed Parameter: Bilateral Stds for the Channel-intensities parsed: ";
			for (int mod_i=0; mod_i<numberOfModalities; mod_i++) {
				parametersOfRun.BilateralModsStds[mod_i] = atof(vectorWithArgs[arg_i + 1 + mod_i].c_str());
				std::cout  << parametersOfRun.BilateralModsStds[mod_i] << " ";
			}
			std::cout << "\n";
			arg_i += numberOfModalities;
		}
		else if (vectorWithArgs[arg_i].compare("-numberOfIterations") == 0) {
			// Number of max iterations to apply the CRF.
			parametersOfRun.MaxIterations = atoi(vectorWithArgs[arg_i+1].c_str());
			arg_i+=1;
			std::cout << "Parsed Parameter: The max number of iterations that the CRF will be applied was parsed to be: " << parametersOfRun.MaxIterations << "\n";
		}
		else {
			std::cout<<"ERROR: Problem when parsing the dense_inference input arg! Unknown option given! Exiting!\n";
			exit(1);
		}
		
	}			

}

// Reads the config file and makes a vector of strings with all the lines of it. Just so that the very same function addPairwiseBilateralToCrfForCorrespondingNumberOfModalities() can process it same way as if read from command line.
void readParamsFromConfigFile(char * configFilePathname, vectorOfStringsType& vectorOfStringArguments) {
    
  std::ifstream infile(configFilePathname);
  std::string newLine;
  while (std::getline(infile, newLine)) {
    vectorOfStringArguments.push_back(newLine);
  }
  infile.close();
}

// Reads the command line arguments and makes a vector of strings with all arguments but the app's name. Just so that the very same function addPairwiseBilateralToCrfForCorrespondingNumberOfModalities() can process it same way as if read from config file.
void readArgvIntoVectorOfStrings(int argc, char* argv[], vectorOfStringsType& vectorOfStringArguments) {
	for (int arg_i=1; arg_i < argc ; arg_i++) {
		vectorOfStringArguments.push_back( std::string(argv[arg_i]) );
		//std::cout<<" DEBUG: COMMAND LINE ARGUMENT: " << vectorOfStringArguments.back() << "\n";
	}
}




int main( int argc, char* argv[]){
	
	ParametersOfRun parametersOfRun;
	
	vectorOfStringsType vectorWithAllModalitiesPaths;
    vectorOfStringsType vectorWithAllProbabilityMapsPaths;
    std::string stringWithTheFolderForThePatientToSaveTheResultsIn = "placeholder";
    std::string stringWithThePrefixOfTheFilenameToSaveResultingSegmMap = "denseCrf3dOutputSegm"; //default name to save resulting segm map
    std::string stringWithThePrefixOfTheFilenameToSaveResultingProbMaps = "denseCrf3dProbMapClass"; //default prefix to save the resulting class prob-maps.
	int numberOfModalities = 0;
	int numberOfForegroundClasses = 0;
	int sizeOfImages[3] = {-1,-1,-1};
	// minIntensity, maxIntensity: All the channels (modalities) need to be at the same intensity range. These variables should be given by the user the corresponding min/max values of the range in their images.
	// Values below or above these limits will be set to the boundary values. 
	// -3 and +3 are the values that were used in our work, Kamnitsas et al, "Multi-Scale 3D Convolutional Neural Networks for Lesion Segmentation in Brain MRI", winning submission to ISLES'15,
	// ...where the intensities of each channel were normalised to 0-mean, 1-std.
	float minIntensity = -3, maxIntensity = +3; 
	
	vectorOfStringsType vectorOfStringArguments;

	if (std::strcmp(argv[1], "-h") == 0) {
		std::cout << "Please check the README that comes with this software for usage instructions.\n"; //printHelpInformation();
		exit(0);
	} else if (std::strcmp(argv[1], "-c") == 0) {
		readParamsFromConfigFile(argv[2], vectorOfStringArguments);
	} else {
		readArgvIntoVectorOfStrings(argc, argv, vectorOfStringArguments);
	}
	
	parseArgs(
			vectorOfStringArguments,
			vectorWithAllModalitiesPaths,
			vectorWithAllProbabilityMapsPaths,
			stringWithTheFolderForThePatientToSaveTheResultsIn,
			stringWithThePrefixOfTheFilenameToSaveResultingSegmMap,
			stringWithThePrefixOfTheFilenameToSaveResultingProbMaps,
			numberOfModalities,
			numberOfForegroundClasses,
			sizeOfImages,
			minIntensity,
			maxIntensity,
			parametersOfRun);
			
	std::cout << "********** Parameters that will be used for this run: **********\n";
	
	std::cout<<"Parameters: numberOfModalities="<<numberOfModalities<<", numberOfForegroundClasses="<<numberOfForegroundClasses<<"\n";
	std::cout<<"Parameters: sizeOfImages= ["<<sizeOfImages[0]<<", "<<sizeOfImages[1]<<", "<<sizeOfImages[2]<<"]\n";
	std::cout<<"Parameters: min Intensity value in the channels=" << minIntensity << ", max intensity value in the channels=" << maxIntensity <<
		" (NOTE: All channels should have been already normalised to the same intensity range)\n";

	parametersOfRun.printParameters();
		
		
	std::cout<<"********** Going to read the modality niftis... **********\n";
	
	InputImageType::Pointer imgOneOfModsItkTypeForReferenceOfResultingImage;
	ImageIOType::Pointer niftiImageIOForReferenceOfResultingImage;// = ImageIOType::New();
    unsigned char * arrayWithMyImWithAllModalities = readNiiModalitiesOfPatientToNeededFormat( vectorWithAllModalitiesPaths, 
																								numberOfModalities,
																								sizeOfImages,
																								minIntensity,
																								maxIntensity,
																								imgOneOfModsItkTypeForReferenceOfResultingImage,
																								niftiImageIOForReferenceOfResultingImage );
    std::cout << "...Done with Image modalities for patient..!\n";

    std::cout << "********** Loading the probability maps and creating the unary potential matrix... **********\n";
    MatrixXf myUnary = readNiiProbMapsOfPatientAndMakeUnaryEnergyMatrix(vectorWithAllProbabilityMapsPaths, numberOfForegroundClasses, sizeOfImages);

    DenseCRF3D crf3d(sizeOfImages[0], sizeOfImages[1], sizeOfImages[2], numberOfForegroundClasses+1);


    std::cout << "********** Setting up Unary and Pairwise potentials, getting ready for inference. **********\n";
    crf3d.setUnaryEnergy( myUnary );

    crf3d.addPairwiseGaussian( parametersOfRun.PosRStd , parametersOfRun.PosCStd, parametersOfRun.PosZStd, new PottsCompatibility( parametersOfRun.PosW ) );

	addPairwiseBilateralToCrfForCorrespondingNumberOfModalities( crf3d, parametersOfRun, arrayWithMyImWithAllModalities);
	
									

    std::cout << "++++++++++++++++++++++++++ Performing Inference ++++++++++++++++++++++++++\n";
    
    MatrixXf probMapsMatrix = crf3d.inference(parametersOfRun.MaxIterations);
    //Save the results
    makeNiiProbMapFromMatrixXfOfProbsAndSave(stringWithTheFolderForThePatientToSaveTheResultsIn,
     										stringWithThePrefixOfTheFilenameToSaveResultingProbMaps,
											probMapsMatrix,
											numberOfForegroundClasses,
											sizeOfImages,
											imgOneOfModsItkTypeForReferenceOfResultingImage,
											niftiImageIOForReferenceOfResultingImage);
											
	VectorXs segmentationVector = crf3d.currentMap(probMapsMatrix);
	
	makeNiiSegmMapFromVectorXsOfLabelsAndSave(stringWithTheFolderForThePatientToSaveTheResultsIn,
											stringWithThePrefixOfTheFilenameToSaveResultingSegmMap,
											segmentationVector, //first row is the background class.
											sizeOfImages, 
											imgOneOfModsItkTypeForReferenceOfResultingImage,
											niftiImageIOForReferenceOfResultingImage);
	
    delete[] arrayWithMyImWithAllModalities;
    std::cout << "++++++++++++++++++++++++++ Done. ++++++++++++++++++++++++++\n";
    
}
