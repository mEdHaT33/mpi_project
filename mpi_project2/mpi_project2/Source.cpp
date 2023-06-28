#include <iostream>
#include <math.h>
#include <stdlib.h>
#include<string.h>
#include<mpi.h>
#include <stdio.h>
#include<sstream>
#include<msclr\marshal_cppstd.h>
#include <ctime>// include this header 
#pragma once
#include <random>
#using <mscorlib.dll>
#using <System.dll>
#using <System.Drawing.dll>
#using <System.Windows.Forms.dll>
using namespace std;
using namespace msclr::interop;

int* inputImage(int* w, int* h, System::String^ imagePath) //put the size of image in w & h
{
	int* input;


	int OriginalImageWidth, OriginalImageHeight;

	//*********************************************************Read Image and save it to local arrayss*************************	
	//Read Image and save it to local arrayss

	System::Drawing::Bitmap BM(imagePath);

	OriginalImageWidth = BM.Width;
	OriginalImageHeight = BM.Height;
	*w = BM.Width;
	*h = BM.Height;
	int *Red = new int[BM.Height * BM.Width];
	int *Green = new int[BM.Height * BM.Width];
	int *Blue = new int[BM.Height * BM.Width];
	input = new int[BM.Height*BM.Width];
	for (int i = 0; i < BM.Height; i++)
	{
		for (int j = 0; j < BM.Width; j++)
		{
			System::Drawing::Color c = BM.GetPixel(j, i);

			Red[i * BM.Width + j] = c.R;
			Blue[i * BM.Width + j] = c.B;
			Green[i * BM.Width + j] = c.G;

			input[i*BM.Width + j] = ((c.R + c.B + c.G) / 3); //gray scale value equals the average of RGB values

		}

	}
	delete[] Red;
	delete[] Green;
	delete[] Blue;

	return input;
}


void createImage(int* image, int width, int height, int index)
{
	System::Drawing::Bitmap MyNewImage(width, height);


	for (int i = 0; i < MyNewImage.Height; i++)
	{
		for (int j = 0; j < MyNewImage.Width; j++)
		{
			//i * OriginalImageWidth + j
			if (image[i*width + j] < 0)
			{
				image[i*width + j] = 0;
			}
			if (image[i*width + j] > 255)
			{
				image[i*width + j] = 255;
			}
			System::Drawing::Color c = System::Drawing::Color::FromArgb(image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j], image[i*MyNewImage.Width + j]);
			MyNewImage.SetPixel(j, i, c);
		}
	}
	MyNewImage.Save("C://Users//Vip B//Desktop//mpi_project2//New folder//OUT//" + index + ".png");
	cout << "result Image Saved " << index << endl;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void parallel_filter(int* imageData, int ImageWidth, int ImageHeight, int rank, int size,int indx, int kernal_size) {


	int total;
	int* temppp = new int[(ImageHeight / size) * ImageWidth];
	int* local_temp = new int[(ImageHeight / size) * ImageWidth];

	int* output = new int[ImageHeight * ImageWidth];
	MPI_Scatter(imageData, (ImageHeight / size) * ImageWidth, MPI_INT, temppp, (ImageHeight / size) * ImageWidth, MPI_INT, 0, MPI_COMM_WORLD);
	int sum=0;
	int ker = (kernal_size - 1) / 2;
	for (int i = 0; i < ImageHeight / size; i++) {
		for (int j = 0; j < ImageWidth; j++) {
			int total = 0;
			int count = 0;

			for (int m = -ker; m <= ker; ++m) {
				for (int n = -ker; n <= ker; ++n) {
					int neighborRow = i + m;
					int neighborCol = j + n;

					// Check if the neighbor pixel is within the image boundaries
					if (neighborRow >= 0 && neighborRow < ImageHeight && neighborCol >= 0 && neighborCol < ImageWidth) {
						total += temppp[neighborRow * ImageWidth + neighborCol];
						count++;
					}
				}
			}

			int sum = total / count;
			local_temp[i * ImageWidth + j] = sum;
		}
	}


	MPI_Gather(local_temp, (ImageHeight / size) * ImageWidth, MPI_INT, output, (ImageHeight / size) * ImageWidth, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {

		createImage(output, ImageWidth, ImageHeight, indx);
	}
	delete[] output;
	delete[] temppp;
	
}


void seq_filter(int* imageData, int ImageWidth, int ImageHeight,int indx, int kernal_size) {
	int ker = (kernal_size - 1) / 2;
	for (int i = 0; i < ImageHeight; i++)
	{
		for (int j = 0; j < ImageWidth; j++)
		{
			// Apply low-pass filter with kernel 1/9 to image[i * width + j]
			int sum = 0;
			for (int m = i - ker; m <= i + ker; m++)
			{
				for (int n = j - ker; n <= j + ker; n++)
				{
					if (m >= 0 && m < ImageHeight && n >= 0 && n < ImageWidth)
					{
						sum += imageData[m * ImageWidth + n];
					}
				}
			}
			imageData[i * ImageWidth + j] = sum / (kernal_size * kernal_size); // Update the pixel value with the filtered value
		}
	}
	createImage(imageData, ImageWidth, ImageHeight, indx);
}


int main(int argc, char** argv) {
    MPI_Init(NULL, NULL);

    int rank, size;
	int kernal_size = 3;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	string img;
    int ImageWidth = 4, ImageHeight = 4;
	for (int i = 1; i <= 4; i++) {
		ostringstream ss;
		ss << "C://Users//Vip B//Desktop//mpi_project2//New folder//" << i << ".png";
		string img = ss.str();

		System::String^ imagePath = marshal_as<System::String^>(img);
		int* imageData = inputImage(&ImageWidth, &ImageHeight, imagePath);

		clock_t start_s, stop_s;
		int totaltime = 0;
		start_s = clock();

		cout << "parallel one" << endl;

		parallel_filter(imageData, ImageWidth, ImageHeight, rank, size,i,kernal_size);

		cout << "Seq one" << endl;

		//seq_filter(imageData, ImageWidth, ImageHeight,i,kernal_size);


		stop_s = clock();
		totaltime += (stop_s - start_s) / double(CLOCKS_PER_SEC) * 1000;
		cout << "time: " << totaltime << endl;


		delete[] imageData;
	
	}
    //std::string img = "C://Users//Vip B//Desktop//mpi_project2//data//16.png";
	MPI_Finalize();
    return 0;
}

