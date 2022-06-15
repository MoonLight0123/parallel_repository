#include <stdio.h>

/*
 * Refactor firstParallel so that it can run on the GPU.
 */

__global__ void firstParallel()
{
  printf("This should be running in parallel.\n");
}

int main()
{
  /*
   * Refactor this call to firstParallel to execute in parallel
   * on the GPU.
   */

  firstParallel<<<3,3>>>();
cudaDeviceSynchronize();
  /*
   * Some code is needed below so that the CPU will wait
   * for the GPU kernels to complete before proceeding.
   */

}

#include <stdio.h>

__global__ void printSuccessForCorrectExecutionConfiguration()
{

  if(threadIdx.x == 1023 && blockIdx.x == 255)
  {
    printf("Success!\n");
  }
}

int main()
{
  /*
   * This is one possible execution context that will make
   * the kernel launch print its success message.
   */

  printSuccessForCorrectExecutionConfiguration<<<256, 1024>>>();

  /*
   * Don't forget kernel execution is asynchronous and you must
   * sync on its completion.
   */

  cudaDeviceSynchronize();
}