/**
 * @author Yunpeng Men
 * @email ypmen@pku.edu.cn
 * @create date 2020-08-11 11:38:57
 * @modify date 2020-08-11 11:38:57
 * @desc [description]
 */

#include "databuffer.h"

#include <fstream>
#include <string.h>
#include <fftw3.h>

using namespace std;

/* (nsamples, nchans) */

template <typename T>
DataBuffer<T>::DataBuffer()
{
	equalized = false;
	isbusy = false;
	closable = false;
	counter = 0;
	nsamples = 0;
	tsamp = 0.;
	nchans = 0;
	npol = 0;
}

template <typename T>
DataBuffer<T>::DataBuffer(const DataBuffer<T> &databuffer)
{
	equalized = databuffer.equalized;
	isbusy = databuffer.isbusy;
	closable = databuffer.closable;
	counter = databuffer.counter;
	nsamples = databuffer.nsamples;
	tsamp = databuffer.tsamp;
	nchans = databuffer.nchans;
	npol = databuffer.npol;
	frequencies = databuffer.frequencies;
	buffer = databuffer.buffer;
}

template <typename T>
DataBuffer<T> & DataBuffer<T>::operator=(const DataBuffer<T> &databuffer)
{
	equalized = databuffer.equalized;
	isbusy = databuffer.isbusy;
	closable = databuffer.closable;
	counter = databuffer.counter;
	nsamples = databuffer.nsamples;
	tsamp = databuffer.tsamp;
	nchans = databuffer.nchans;
	npol = databuffer.npol;
	frequencies = databuffer.frequencies;
	buffer = databuffer.buffer;

	return *this;    
}

template <typename T>
DataBuffer<T>::DataBuffer(long int ns, int nc)
{
	counter = 0;
	resize(ns, nc, 1);
	tsamp = 0.;

	equalized = false;
	isbusy = false;
	closable = false;
}

template <typename T>
DataBuffer<T>::DataBuffer(long int ns, int nc, int np)
{
 	counter = 0;
        resize(ns, nc, np);
        tsamp = 0.;

        equalized = false;
        isbusy = false;
        closable = false;
}

template <typename T>
DataBuffer<T>::~DataBuffer(){}

template <typename T>
void DataBuffer<T>::prepare(DataBuffer<T> &databuffer)
{
	equalized = databuffer.equalized;
	nsamples = databuffer.nsamples;
	nchans = databuffer.nchans;
	npol = databuffer.npol;
	
	resize(nsamples, nchans, npol);

	tsamp = databuffer.tsamp;
	frequencies = databuffer.frequencies;
}

template <typename T>
DataBuffer<T> * DataBuffer<T>::get_pol(DataBuffer<T> &databuffer, int ipol)
{
        equalized = databuffer.equalized;
        nsamples = databuffer.nsamples;
	nchans = databuffer.nchans;
        npol = 1;

        resize(nsamples, nchans, npol);

	tsamp = databuffer.tsamp;
        frequencies = databuffer.frequencies;

        for (long int i=0; i<nsamples; i++)
        {
	        for (long int j=0; j<nchans; j++)
                {
			buffer[i*nchans+j] = databuffer.buffer[i*databuffer.npol*databuffer.nchans+ipol*databuffer.nchans+j];
                }

        }
	return this;
}

template <typename T>
DataBuffer<T> * DataBuffer<T>::four2one(DataBuffer<T> &databufferI, DataBuffer<T> &databufferQ, DataBuffer<T> &databufferU, DataBuffer<T> &databufferV)
{
  //buffer = databuffer.buffer;
  equalized = databufferI.equalized;
  nsamples = databufferI.nsamples;
  nchans = databufferI.nchans;
  npol = 4;
  counter = databufferI.counter;
  equalized = databufferI.equalized;
  isbusy = databufferI.isbusy;
  resize(nsamples, nchans, npol);

  tsamp = databufferI.tsamp;
  frequencies = databufferI.frequencies;
  
	for (long int i=0; i<nsamples; i++)
	  {
	    copy(databufferI.buffer.begin() + i*nchans, databufferI.buffer.begin() + (i+1)*nchans+1, buffer.begin()+i*npol*nchans);
	    copy(databufferQ.buffer.begin() + i*nchans, databufferQ.buffer.begin() + (i+1)*nchans+1, buffer.begin()+i*npol*nchans + 1*nchans);
	    copy(databufferU.buffer.begin() + i*nchans, databufferU.buffer.begin() + (i+1)*nchans+1, buffer.begin()+i*npol*nchans + 2*nchans);
	    copy(databufferV.buffer.begin() + i*nchans, databufferV.buffer.begin() + (i+1)*nchans+1, buffer.begin()+i*npol*nchans + 3*nchans);
	  }
  
        return this;
};

template <typename T>
DataBuffer<T> * DataBuffer<T>::run(DataBuffer<T> &databuffer)
{
	buffer = databuffer.buffer;

	counter += nsamples;

	databuffer.isbusy = false;
	isbusy = true;

	return this;
};

template <typename T>
DataBuffer<T> * DataBuffer<T>::filter(DataBuffer<T> &databuffer)
{
	counter += nsamples;

	databuffer.isbusy = true;

	return databuffer.get();
};

template <typename T>
void DataBuffer<T>::open()
{
	buffer.clear();
	buffer.resize(nsamples*nchans*npol, 0.);
}

template <typename T>
void DataBuffer<T>::close()
{
	buffer.clear();
	buffer.shrink_to_fit();
}

template <typename T>
void DataBuffer<T>::dump2txt(const string fname)
{
  string st = "-Pol";
	ofstream outfile;
	for (int k=0; k<npol; k++)
	  {
	    outfile.open(fname+st+to_string(k), ofstream::out);

	    for (long int i=0; i<nsamples; i++)
	      {
		for (long int j=0; j<nchans; j++)
		  {
		    outfile<<buffer[i*npol*nchans+k*nchans+j]<<" ";
		  }
		outfile<<endl;
	      }

	    outfile.close();
	  }
}

template <typename T>
void DataBuffer<T>::dump2bin(const string fname)
{
	ofstream outfile;
	outfile.open(fname, ofstream::binary);

	outfile.write((char *)(&buffer[0]), sizeof(T)*nsamples*nchans);

	outfile.close();
}

template <typename T>
void DataBuffer<T>::dump(const string fname)
{
	ofstream outfile;
	outfile.open(fname, ios::binary|ios::app);

	outfile.write((char *)(&buffer[0]), sizeof(T)*nsamples*nchans);

	outfile.close();
}
//template <typename T>
//void DataBuffer<T>::resize(long int ns, int nc)
//{
//	nsamples = ns;
//	nchans = nc;
//	buffer.resize(nsamples*nchans*npol, 0.);
//	frequencies.resize(nchans, 0.);
//}

template <typename T>
void DataBuffer<T>::resize(long int ns, int nc, int np)
{
        nsamples = ns;
        nchans = nc;
	npol = np;
        buffer.resize(nsamples*nchans*npol, 0.);
        frequencies.resize(nchans, 0.);
}

template <typename T>
void DataBuffer<T>::get_mean_rms(vector<T> &mean, vector<T> &var)
{
	mean.resize(nchans, 0.);
	var.resize(nchans, 0.);
	for (long int i=0; i<nsamples; i++)
	{
	  for (int k=0; k<npol; k++)
	    {
		for (long int j=0; j<nchans; j++)
		{
			mean[j] += buffer[i*npol*nchans+k*nchans+j];
			var[j] += buffer[i*npol*nchans+k*nchans+j]*buffer[i*npol*nchans+k*nchans+j];
		}
	    }
	}

	for (long int j=0; j<nchans; j++)
	{
		mean[j] /= nsamples;
		var[j] /= nsamples;
		var[j] -= mean[j]*mean[j];
	}
}

template class DataBuffer<char>;
template class DataBuffer<unsigned char>;
template class DataBuffer<short>;
template class DataBuffer<float>;
template class DataBuffer<double>;
template class DataBuffer<complex<float>>;
template class DataBuffer<complex<double>>;
