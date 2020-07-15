/*
 * Copyright 2016 <Admobilize>
 * MATRIX Labs  [http://creator.matrix.one]
 * This file is part of MATRIX Creator HAL
 *
 * MATRIX Creator HAL is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cmath>
#include <map>
#include <string>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <complex>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/FFT>

#include "cpp/driver/cross_correlation.h"
#include "cpp/driver/microphone_array.h"
#include "cpp/driver/direction_of_arrival_music.h"
#include "cpp/driver/microphone_array_location.h"

#define PI 3.14159265358979323846
using namespace std;
using namespace Eigen;

//typedef Matrix<float, 8, 128> Matrix8by128f;
//typedef Matrix<float, 128, 8> Matrix128by8f;
//typedef Matrix<float, 128, 128> Matrix128by128f;
typedef Matrix<float, 8, 8> Matrix8by8f;

typedef Matrix<complex<float>, 8, 8> MatrixC8by8f;
typedef Matrix<complex<float>, 1, 8> MatrixC1by8f;
//typedef Matrix<MatrixC8by8f, 128, 1> MatrixCC128by1f;

typedef Matrix<float, 8, 7> Matrix8by7f;
typedef Matrix<float, 7, 8> Matrix7by8f;

typedef Matrix<complex<float>, 8, 7> MatrixC8by7f;
typedef Matrix<complex<float>, 7, 8> MatrixC7by8f;
//typedef Matrix<float, 128, 90> Matrix128by90f;
//typedef Matrix<float, 128, 38> Matrix128by38f;
typedef Matrix<float, 8, 1> Matrix8by1f;
typedef Matrix<complex<float>, 8, 1> MatrixC8by1f;
typedef Matrix<complex<float>, 8, 5> MatrixC8by5f;
typedef Matrix<float, 1, 1> Matrix1by1f;
typedef Matrix<complex<float>, 1, 1> MatrixC1by1f;
typedef Matrix<Matrix8by1f, 360, 1> MatrixN360by1f;
typedef Matrix<float, 360, 1> Matrix360by1f;
typedef Matrix<MatrixC8by1f, 360, 1> MatrixC360by1f;
typedef Matrix<MatrixC8by1f, 1, 360> MatrixC1by360f;

MatrixC360by1f steeringVector_;

namespace matrix_hal {

DirectionOfArrivalMusic::DirectionOfArrivalMusic(MicrophoneArray &mics) : mics_(mics) {}

bool DirectionOfArrivalMusic::Init() {
    length_ = mics_.NumberOfSamples();
    corr_ = new CrossCorrelation();
    corr_->Init(mics_.NumberOfSamples());
    current_mag_.resize(4);
    current_index_.resize(4);
    buffer_1D_.resize(mics_.Channels() * mics_.NumberOfSamples());
    buffer_2D_.resize(mics_.Channels());
    mic_direction_ = 0;
    azimutal_angle_ = 0;
    polar_angle_ = 0;
    snap_shots_ = (int)(mics_.SamplingRate()/mics_.NumberOfSamples())+1; // 125.0; // 1000 ms of data
    file_count_ = 0.0;
    for (uint16_t c = 0; c < mics_.Channels(); c++) {
        buffer_2D_[c] = &buffer_1D_[c * mics_.NumberOfSamples()];
    }
    steeringVectorCalculate();
    return true;
}

void DirectionOfArrivalMusic::steeringVectorCalculate(){
    const complex<double> l(0,1);
    const float radius_matrix_voice = 0.0383; // use 0.054 for matrix creator
    for(int i = 0 ; i < 360 ; i ++){
        MatrixC8by1f steeringVec_;
        for(int j = 0 ; j < 8 ; j++){
            steeringVec_(j,0) = std::exp(l*( ((2*PI*mics_.SamplingRate()*radius_matrix_voice)/340.0)*
                                            //cos( (i*PI/180.0- (2*PI*PI*j)/(8.0*180.0)) ) ) );
                                            cos( (i*PI/180.0- (2*PI*j)/(8.0)) ) ) );
        }
        steeringVector_(i,0) = MatrixC8by1f(steeringVec_);
    }
}

int DirectionOfArrivalMusic::getAbsDiff(int index) {
    if (index < length_ / 2) {
        return index;
    }
    return length_ - 1 - index;
}

void DirectionOfArrivalMusic::Calculate() {
    //Create the Matrix
    const uint16_t numCh = 8;
    const uint32_t numS = 512;
    MatrixC8by8f musicFinalMatrix;
    
    for(int i = 0; i < snap_shots_; i++){
        
        FFT<float> fft;
        Matrix<float, numCh, numS> musicMatrix;
        Matrix<complex<float>, numCh, numS> outMatrix;
        outMatrix.setZero(numCh, numS);
        
        //Putting the raw Data into buffers
        for (uint32_t s = 0; s < mics_.NumberOfSamples(); s++) {
            for (uint16_t c = 0; c < mics_.Channels(); c++) {
                buffer_2D_[c][s] = mics_.Raw(s, c); //Raw data (not delayed data)
                musicMatrix(c,s) = buffer_2D_[c][s];
            }
        }
        
        //FFT Operations
        for (int k = 0; k < musicMatrix.rows(); k++) {
            Matrix<complex<float>, 1, numS> tmpOut;
            fft.fwd(tmpOut, musicMatrix.row(k));
            outMatrix.row(k) = tmpOut;
        }
        
        for (int k = 0; k < musicMatrix.cols(); k++) {
            Matrix<complex<float>, 8, 1> tmpOut;
            fft.fwd(tmpOut, outMatrix.col(k));
            outMatrix.col(k) = tmpOut;
        }
        
        //Reading again Mic Values
        mics_.Read();
        
        //Calculating the average matrix values
        if(i == 0 ){
            Matrix<complex<float>, numS, numCh> outMatrixAdj(outMatrix.adjoint());
            musicFinalMatrix = outMatrix*(outMatrixAdj);
        }else{
            MatrixC8by8f tmpMtx;
            Matrix<complex<float>, numS, numCh> outMatrixAdj(outMatrix.adjoint());
            tmpMtx = musicFinalMatrix + (outMatrix*(outMatrixAdj));
            musicFinalMatrix = tmpMtx;
        }
        //cout << "----------------------" << endl << musicFinalMatrix << endl << "----------------------" << endl;
    } //end of for loop
    
    //musicFinalMatrix.setZer(8,8);
    musicFinalMatrix = (1.0/(1.0*snap_shots_))*musicFinalMatrix;
    
    //For Complex Mtarix
    ComplexEigenSolver<MatrixC8by8f> es;
    es.compute(musicFinalMatrix);
    
    cout << "The eigenvalues of Correlation Matrix are:" << endl << es.eigenvalues() << endl;
    
    MatrixC8by7f Qmat;
    
    for(int l=0; l< 7 ; l++)
        Qmat.col(l) = es.eigenvectors().col(l);
    
    Matrix360by1f phaseSpectrum_;
    
    //ofstream myfile;
    
    //myfile.open ("example_" + std::to_string(file_count_) + ".dat");
    MatrixC7by8f QmatAdj(Qmat.adjoint());
    
    for(int i = 0; i < 360; i++){
        MatrixC1by8f stVecAdj(steeringVector_(i,0).adjoint());
        MatrixC1by1f val = stVecAdj*Qmat*QmatAdj*steeringVector_(i,0);
        phaseSpectrum_(i,0) = 1.0/(val.norm());
        //myfile << i << " " << 1.0/(val.norm()) << endl;
    }
    
    //myfile.close();
    file_count_ += 1;
    
    //get location of maximum
    MatrixXf::Index maxRow, maxCol;
    float max = phaseSpectrum_.maxCoeff(&maxRow, &maxCol);
    //cout << phaseSpectrum_ << endl;
    
    cout << "Value :" << max << "; Max Degree : " << maxRow << endl;
    
    if ( max > 0.05) {
        mic_direction_ = maxRow/45;
        azimutal_angle_ = atan2(micarray_location_voice[mic_direction_][1],
                                micarray_location_voice[mic_direction_][0]);
        //polar_angle_ = fabs(index) * M_PI / 2.0 / float(max_tof - 1);
    }
}

};  // namespace matrix_hal
