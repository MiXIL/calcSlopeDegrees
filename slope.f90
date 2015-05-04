c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c Slope calculation routine. Designed to be called from python
c
c Daniel Clewley (daniel.clewley@gmail.com) - 27/06/2013
c
c Copyright 2014 Daniel Clewley & Jane Whitcomb.
c 
c Permission is hereby granted, free of charge, to any person 
c obtaining a copy of this software and associated documentation 
c files (the "Software"), to deal in the Software without restriction, 
c including without limitation the rights to use, copy, modify, 
c merge, publish, distribute, sublicense, and/or sell copies of the 
c Software, and to permit persons to whom the Software is furnished 
c to do so, subject to the following conditions:
c 
c The above copyright notice and this permission notice shall be 
c included in all copies or substantial portions of the Software.
c 
c THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
c EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
c OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
c IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
c ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF 
c CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
c WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      subroutine slope(inBlock, xSizes, ySizes, width, height, zScale
     1    , outBlock)
c     slope routine - standard slope calculation

c     f2py will figure out width and height are the array size and set
c     automatically so these don't need to be passed in.

      implicit none

      integer x, y
      real pi
      real dx, dy, dzx, dzy
      real nx, ny, nz
      real slopeRad, slopeDeg
      integer width, height
      real zScale

c     Array of x and y sizes is 2 pixels smaller than data array (for overlap)      
      real inBlock(height,width)
      real outBlock(height,width)
      real xSizes(height, width)
      real ySizes(height, width)

      pi = 3.14159265

c     Tell f2py that outBlock is the output

Cf2py intent(out) outBlock

      do x=2,width-1
        do y=2,height-1
c           Get pixel size
            dx = 2 * xSizes(y,x)
            dy = 2 * ySizes(y,x)
c            Calculate difference in elevation
            dzx = (inBlock(y,x-1) - inBlock(y,x+1))*zScale
            dzy = (inBlock(y-1,x) - inBlock(y+1,x))*zScale

c            Find normal vector to the plane
            nx = -1 * dy * dzx
            ny = -1 * dx * dzy
            nz = dx * dy

            slopeRad = acos(nz / sqrt(nx**2 + ny**2 + nz**2))
            slopeDeg = (180. / pi) * slopeRad

            outBlock(y,x) = slopeDeg
        end do

      end do

      return 
      end

      subroutine slopeplane(inBlock, xSizes, ySizes, width, height,
     1       zScale, outBlock, winSize)
c     slopeplane routine - slope calculation using least squares
c     to fit a plane to elevation.

c     f2py will figure out width and height are the array size and set
c     automatically so these don't need to be passed in.

      implicit none

      external sgels
      intrinsic int, min

      integer x, y, i, j, indx
      integer winSize, winOffset
      real pi
      real dx, dy, dzx, dzy
      real nx, ny, nz
      real slopeRad, slopeDeg
      integer nRows, nCols
      integer width, height
      real zScale
      integer lwork, info

c     Array of x and y sizes is 2 pixels smaller than data array (for overlap)      
      real inBlock(height,width)
      real outBlock(height,width)
      real xSizes(height, width)
      real ySizes(height, width)
      parameter (nCols = 3)
      real AMat(winSize**2,nCols)
      real zVec(winSize**2)
      real work(winSize**2*nCols*2)
      nRows = winSize**2

      winOffset=winSize/2

      pi = 3.14159265

c     Tell f2py that outBlock is the output

Cf2py intent(out) outBlock

      indx=1
      do x=winOffset,width-winOffset
        do y=winOffset,height-winOffset
c           Get pixel size
            dx = winSize * xSizes(y,x)
            dy = winSize * ySizes(y,x)

c           Calculate difference in elevation

c           Set up matrices for plane fitting

            indx = 1
            do i=(-1*winOffset),(winOffset)
                do j=(-1*winOffset),(winOffset)

                    AMat(indx,1) = (i*xSizes(y,x))
                    AMat(indx,2) = (j*ySizes(y,x))
                    AMat(indx,3) = 1

                    zVec(indx) = inBlock(y+j,x+i)*zScale

                    indx=indx+1

                end do
            end do

c           Fit plane
            dzx=1.
            dzy=1.

            lwork=-1
c           Get optimal value for work array
            call sgels('N', nRows,nCols,1,AMat,nRows,zVec,nRows,work,
     1         lwork, info)
            lwork=int(work(1))
c           Actual fitting
            call sgels('N', nRows,nCols,1,AMat,nRows,zVec,nRows,work,
     1         lwork, info)


c           Calculate dzx and dzy from coefficients of fit
            dzx = zVec(1) * dx
            dzy = zVec(2) * dy

c            Find normal vector to the plane
            nx = -1 * dy * dzx
            ny = -1 * dx * dzy
            nz = dx * dy

            slopeRad = acos(nz / sqrt(nx**2 + ny**2 + nz**2))
            slopeDeg = (180. / pi) * slopeRad

            outBlock(y,x) = slopeDeg
        end do

      end do

      return 
      end


 
