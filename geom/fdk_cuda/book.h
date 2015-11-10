/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * NVIDIA Corporation and its licensors retain all intellectual property and
 * proprietary rights in and to this software and related documentation.
 * Any use, reproduction, disclosure, or distribution of this software
 * and related documentation without an express license agreement from
 * NVIDIA Corporation is strictly prohibited.
 *
 * Please refer to the applicable NVIDIA end user license agreement (EULA)
 * associated with this source code for terms and conditions that govern
 * your use of this NVIDIA software.
 *
 */


#ifndef __BOOK_H__
#define __BOOK_H__
#include <stdio.h>

static void HandleError( cudaError_t err,
        const char *file,
        int line ) {
    if (err != cudaSuccess) {
        mexPrintf( "%s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))


#define HANDLE_NULL( a ) {if (a == NULL) { \
                            mexPrintf( "Host memory failed in %s at line %d\n", \
                                    __FILE__, __LINE__ ); \
                                    exit( EXIT_FAILURE );}}

template< typename T >
        void swap( T& a, T& b ) {
    T t = a;
    a = b;
    b = t;
}


void* big_random_block( int size ) {
    unsigned char *data = (unsigned char*)malloc( size );
    HANDLE_NULL( data );
    for (int i=0; i<size; i++)
        data[i] = rand();
    
    return data;
}

int* big_random_block_int( int size ) {
    int *data = (int*)malloc( size * sizeof(int) );
    HANDLE_NULL( data );
    for (int i=0; i<size; i++)
        data[i] = rand();
    
    return data;
}


// a place for common kernels - starts here

__device__ unsigned char value( float n1, float n2, int hue ) {
    if (hue > 360)      hue -= 360;
    else if (hue < 0)   hue += 360;
    
    if (hue < 60)
        return (unsigned char)(255 * (n1 + (n2-n1)*hue/60));
    if (hue < 180)
        return (unsigned char)(255 * n2);
    if (hue < 240)
        return (unsigned char)(255 * (n1 + (n2-n1)*(240-hue)/60));
    return (unsigned char)(255 * n1);
}

__global__ void float_to_color( unsigned char *optr,
        const float *outSrc ) {
    // map from threadIdx/BlockIdx to pixel position
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = x + y * blockDim.x * gridDim.x;
    
    float l = outSrc[offset];
    float s = 1;
    int h = (180 + (int)(360.0f * outSrc[offset])) % 360;
    float m1, m2;
    
    if (l <= 0.5f)
        m2 = l * (1 + s);
    else
        m2 = l + s - l * s;
    m1 = 2 * l - m2;
    
    optr[offset*4 + 0] = value( m1, m2, h+120 );
    optr[offset*4 + 1] = value( m1, m2, h );
    optr[offset*4 + 2] = value( m1, m2, h -120 );
    optr[offset*4 + 3] = 255;
}

__global__ void float_to_color( uchar4 *optr,
        const float *outSrc ) {
    // map from threadIdx/BlockIdx to pixel position
    int x = threadIdx.x + blockIdx.x * blockDim.x;
    int y = threadIdx.y + blockIdx.y * blockDim.y;
    int offset = x + y * blockDim.x * gridDim.x;
    
    float l = outSrc[offset];
    float s = 1;
    int h = (180 + (int)(360.0f * outSrc[offset])) % 360;
    float m1, m2;
    
    if (l <= 0.5f)
        m2 = l * (1 + s);
    else
        m2 = l + s - l * s;
    m1 = 2 * l - m2;
    
    optr[offset].x = value( m1, m2, h+120 );
    optr[offset].y = value( m1, m2, h );
    optr[offset].z = value( m1, m2, h -120 );
    optr[offset].w = 255;
}


#if _WIN32
//Windows threads.
#include <windows.h>

typedef HANDLE CUTThread;
typedef unsigned (WINAPI *CUT_THREADROUTINE)(void *);

#define CUT_THREADPROC unsigned WINAPI
#define  CUT_THREADEND return 0

#else
//POSIX threads.
#include <pthread.h>

typedef pthread_t CUTThread;
typedef void *(*CUT_THREADROUTINE)(void *);

#define CUT_THREADPROC void
#define  CUT_THREADEND
#endif

//Create thread.
CUTThread start_thread( CUT_THREADROUTINE, void *data );

//Wait for thread to finish.
void end_thread( CUTThread thread );

//Destroy thread.
void destroy_thread( CUTThread thread );

//Wait for multiple threads.
void wait_for_threads( const CUTThread *threads, int num );

#if _WIN32
//Create thread
CUTThread start_thread(CUT_THREADROUTINE func, void *data){
    return CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)func, data, 0, NULL);
}

//Wait for thread to finish
void end_thread(CUTThread thread){
    WaitForSingleObject(thread, INFINITE);
    CloseHandle(thread);
}

//Destroy thread
void destroy_thread( CUTThread thread ){
    TerminateThread(thread, 0);
    CloseHandle(thread);
}

//Wait for multiple threads
void wait_for_threads(const CUTThread * threads, int num){
    WaitForMultipleObjects(num, threads, true, INFINITE);
    
    for(int i = 0; i < num; i++)
        CloseHandle(threads[i]);
}

#else
//Create thread
CUTThread start_thread(CUT_THREADROUTINE func, void * data){
    pthread_t thread;
    pthread_create(&thread, NULL, func, data);
    return thread;
}

//Wait for thread to finish
void end_thread(CUTThread thread){
    pthread_join(thread, NULL);
}

//Destroy thread
void destroy_thread( CUTThread thread ){
    pthread_cancel(thread);
}

//Wait for multiple threads
void wait_for_threads(const CUTThread * threads, int num){
    for(int i = 0; i < num; i++)
        end_thread( threads[i] );
}

#endif

void enum_gpu( ) {
    cudaDeviceProp  prop;
    
    int count;
    HANDLE_ERROR( cudaGetDeviceCount( &count ) );
    for (int i=0; i< count; i++) {
        HANDLE_ERROR( cudaGetDeviceProperties( &prop, i ) );
        mexPrintf( "   --- General Information for device %d ---\n", i );
        mexPrintf( "Name:  %s\n", prop.name );
        mexPrintf( "Compute capability:  %d.%d\n", prop.major, prop.minor );
        mexPrintf( "Clock rate:  %d\n", prop.clockRate );
        mexPrintf( "Device copy overlap:  " );
        if (prop.deviceOverlap)
            mexPrintf( "Enabled\n" );
        else
            mexPrintf( "Disabled\n");
        mexPrintf( "Kernel execution timeout :  " );
        if (prop.kernelExecTimeoutEnabled)
            mexPrintf( "Enabled\n" );
        else
            mexPrintf( "Disabled\n" );
        
        mexPrintf( "   --- Memory Information for device %d ---\n", i );
        mexPrintf( "Total global mem:  %ld\n", prop.totalGlobalMem );
        mexPrintf( "Total constant Mem:  %ld\n", prop.totalConstMem );
        mexPrintf( "Max mem pitch:  %ld\n", prop.memPitch );
        mexPrintf( "Texture Alignment:  %ld\n", prop.textureAlignment );
        
        mexPrintf( "   --- MP Information for device %d ---\n", i );
        mexPrintf( "Multiprocessor count:  %d\n",
                prop.multiProcessorCount );
        mexPrintf( "Shared mem per mp:  %ld\n", prop.sharedMemPerBlock );
        mexPrintf( "Registers per mp:  %d\n", prop.regsPerBlock );
        mexPrintf( "Threads in warp:  %d\n", prop.warpSize );
        mexPrintf( "Max threads per block:  %d\n",
                prop.maxThreadsPerBlock );
        mexPrintf( "Max thread dimensions:  (%d, %d, %d)\n",
                prop.maxThreadsDim[0], prop.maxThreadsDim[1],
                prop.maxThreadsDim[2] );
        mexPrintf( "Max grid dimensions:  (%d, %d, %d)\n",
                prop.maxGridSize[0], prop.maxGridSize[1],
                prop.maxGridSize[2] );
        mexPrintf( "\n" );
    }
}



#endif  // __BOOK_H__
