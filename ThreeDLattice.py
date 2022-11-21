import numpy as np 
import numpy.linalg
import scipy.linalg


def ThreeDLattice(nx,ny,nz,diagterm,frontterm,upterm,leftterm):


    #Takes size input and unique terms and returns total matrix, as well as eigenvalues and eigenvectors for the stationary solution
    
    
    #idea: d_x will give right and left values -/+1 cell away from diagional on N_x by N_x matrix
    #d_y will give up and down values -/+ N_x entries away from diagional on (N_x*N_y) by (N_x*N_y) matrix
    #d_z will give front and back values -/+ N_x*N_y entries away from diagional on (N_x*N_y*N_z) by (N_x*N_y*N_z) matrix
    
    #Total matrix is intexed as zyx where x is looped over, and then y is looped over, and then z
    
    #for dimensions (n by m by w)
    #first row would look like [000, 00x1, ..., 00x_n, 0(y_1)0, 0(y_1)x_1, ... 0(y_1)x_n, ... 0(y_m)x_n, (z_1)00, (z_1)0x_1, ... (z_1)y_m(x_n), ..., z_w(y_m)x_n] 
    
    #diagterms = 1. 
    #nx = 3
    #ny = 3 
    #nz = 4

    #frontterm = 2. 
    backterm = -frontterm

    #upterm = 5.
    downterm = -downterm

    #leftterm = 7.
    rightterm=-leftterm

    diag = diagterm*np.identity(nx*ny*nz)

    front = frontterm*np.identity(nx*ny*nz)
    front = np.roll(front, nx*ny,axis=1)

    back = backterm*np.identity(nx*ny*nz)
    back = np.roll(back, -nx*ny,axis=1)

    up = upterm*np.identity(nx*ny)
    up = np.roll(up,nx,axis=1)
    #need to have nz copies of 'up' in the following function
    uppertry = scipy.linalg.block_diag(up,up,up,up)

    down = downterm*np.identity(nx*ny)
    down = np.roll(down,-nx,axis=1)
    #need to have nz copies of 'down' in the following function
    downtry = scipy.linalg.block_diag(down,down,down,down)

    left = leftterm*np.identity(nx)
    left = np.roll(left,-1, axis=1)
    #need to have nz*ny copies of 'left' in the following function
    lefttry = scipy.linalg.block_diag(left,left,left,left,left,left,left,left,left,left,left,left)

    right = rightterm*np.identity(nx)
    right = np.roll(right,1,axis=1)
    #need to have nz*ny copies of 'right' in the following function
    righttry = scipy.linalg.block_diag(right,right,right,right,right,right,right,right,right,right,right,right)


    total = diag + front+back+uppertry+downtry+lefttry+righttry


    eigenval,eigenvect = np.linalg.eig(total)
    
    return total,eigenval,eigenvect



