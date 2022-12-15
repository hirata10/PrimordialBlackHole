import numpy as np 
import numpy.linalg
import scipy.linalg


def ThreeDLattice(nx,ny,nz,m,Ex,Ey,Ez):


    #Takes size input, m, and lattice spacing and returns total matrix, as well as eigenvalues and eigenvectors 
    #WITHOUT WILSON LOOP CORRECTION!!
    
    
    
    #idea: d_x will give right and left values -/+4 cell away from diagional on 4N_x by 4N_x matrix
    #d_y will give up and down values -/+ 4N_x entries away from diagional on 4(N_x*N_y) by 4(N_x*N_y) matrix
    #d_z will give front and back values -/+ 4N_x*N_y entries away from diagional on 4(N_x*N_y*N_z) by 4(N_x*N_y*N_z) matrix
    
    #Total matrix is indexed as zyxa where a (dirac index [0123]) is looped over, x is looped over, then y is looped over, and then z
    
    #for dimensions (n by m by w)
    #first row would look like [0000,0001, 0002,0003 00(x_1)0, ..., 00(x_n)3, 0(y_1)00,... 0(y_1)(x_1)0, ... 0(y_1)(x_n)3, ... 0(y_m)x_n, (z_1)00, (z_1)0x_1, ... (z_1)y_m(x_n), ..., z_w(y_m)(x_n)3] 

    
    #dt = complex(omega,0)*np.diag([complex(0,1),complex(0,1),complex(0,-1),complex(0,-1)]) - m*np.identity(4)
    dt = -m*np.identity(4)
    
    #nx = 3
    #ny = 3 
    #nz = 4

    bt =  np.array([[0,0,1.,0],[0,0,0,-1.],[1.,0,0,0],[0,-1.,0,0]])/(2*Ez)
    ft =  np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez)

    ut  = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey)
    dot = np.array([[0,0,0,complex(0,1)],[0,0,complex(0,-1),0],[0,complex(0,1),0,0],[complex(0,-1),0,0,0]])/(2*Ey)
    
    
    lt = np.array([[0,0,0,-1],[0,0,-1,0],[0,-1,0,0],[-1,0,0,0]])/(2*Ex)
    rt = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex)
    
    
    a = nx*ny*nz
    ai = 0 
    diag = np.array(dt)
    front = np.array(ft)
    back = np.array(bt)
    
    beta = np.array([[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]])
    beta_blocks = np.array(beta)
    while (ai<a-1) ==True:
        diag = scipy.linalg.block_diag(diag,dt)
        front = scipy.linalg.block_diag(front,ft)
        back = scipy.linalg.block_diag(back,bt)
        beta_blocks = scipy.linalg.block_diag(beta_blocks, beta)
        ai+=1

    front = np.roll(front, 4*nx*ny,axis=1)
    back = np.roll(back, -4*nx*ny,axis=1)
    
    
    a = nx*ny
    ai = 0 
    up = np.array(ut)
    down = np.array(dot)
    while (ai<a-1) ==True:
        up = scipy.linalg.block_diag(up,ut)
        down = scipy.linalg.block_diag(down,dot)
        ai+=1

    up = np.roll(up,4*nx,axis=1)
    down = np.roll(down,-4*nx,axis=1)
    
    
    a = nz
    ai = 0 
    uppertry = np.array(up)
    downtry = np.array(down)
    while (ai<a-1) ==True:
        uppertry = scipy.linalg.block_diag(uppertry,up)
        downtry = scipy.linalg.block_diag(downtry,down)
        ai+=1

    a = nx
    ai = 0 
    left = np.array(lt)
    right = np.array(rt)
    while (ai<a-1) ==True:
        left = scipy.linalg.block_diag(left,lt)
        right = scipy.linalg.block_diag(right,rt)
        ai+=1

    left = np.roll(left,-4,axis=1)
    right = np.roll(right,4,axis=1)
    
    
    a = ny*nz
    ai = 0 
    lefttry = np.array(left)
    righttry = np.array(right)
    while (ai<a-1) ==True:
        lefttry = scipy.linalg.block_diag(lefttry,left)
        righttry = scipy.linalg.block_diag(righttry,right)
        ai+=1

    
    total = diag+front+back+uppertry+downtry+lefttry+righttry
    
    beta_total = np.dot(beta_blocks,total)
    #checked and this is hermitian which is great so now can use .eigh

    eigenval,eigenvect = np.linalg.eigh(total)
    
    return beta_total,eigenval,eigenvect


def ThreeDLatticeW(nx,ny,nz,m,Ex,Ey,Ez,rx,ry,rz):


    #Takes size input, m, and lattice spacing and returns total matrix, as well as eigenvalues and eigenvectors 
    #WITH WILSON LOOP CORRECTION!!
    
    
    
    #idea: d_x will give right and left values -/+4 cell away from diagional on 4N_x by 4N_x matrix
    #d_y will give up and down values -/+ 4N_x entries away from diagional on 4(N_x*N_y) by 4(N_x*N_y) matrix
    #d_z will give front and back values -/+ 4N_x*N_y entries away from diagional on 4(N_x*N_y*N_z) by 4(N_x*N_y*N_z) matrix
    
    #Total matrix is indexed as zyxa where a (dirac index [0123]) is looped over, x is looped over, then y is looped over, and then z
    
    #for dimensions (n by m by w)
    #first row would look like [0000,0001, 0002,0003 00(x_1)0, ..., 00(x_n)3, 0(y_1)00,... 0(y_1)(x_1)0, ... 0(y_1)(x_n)3, ... 0(y_m)x_n, (z_1)00, (z_1)0x_1, ... (z_1)y_m(x_n), ..., z_w(y_m)(x_n)3] 

    
    #dt = complex(omega,0)*np.diag([complex(0,1),complex(0,1),complex(0,-1),complex(0,-1)]) - m*np.identity(4)
    dt = -m*np.identity(4) + (rx+ry+rz)*np.identity(4)
    
    #nx = 3
    #ny = 3 
    #nz = 4

    bt =  np.array([[0,0,1.,0],[0,0,0,-1.],[1.,0,0,0],[0,-1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
    ft =  np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)

    ut  = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
    dot = np.array([[0,0,0,complex(0,1)],[0,0,complex(0,-1),0],[0,complex(0,1),0,0],[complex(0,-1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
    
    
    lt = np.array([[0,0,0,-1],[0,0,-1,0],[0,-1,0,0],[-1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
    rt = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
    
    
    a = nx*ny*nz
    ai = 0 
    diag = np.array(dt)
    front = np.array(ft)
    back = np.array(bt)
    
    beta = np.array([[-1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]])
    beta_blocks = np.array(beta)
    while (ai<a-1) ==True:
        diag = scipy.linalg.block_diag(diag,dt)
        front = scipy.linalg.block_diag(front,ft)
        back = scipy.linalg.block_diag(back,bt)
        beta_blocks = scipy.linalg.block_diag(beta_blocks, beta)
        ai+=1

    front = np.roll(front, 4*nx*ny,axis=1)
    back = np.roll(back, -4*nx*ny,axis=1)
    
    
    a = nx*ny
    ai = 0 
    up = np.array(ut)
    down = np.array(dot)
    while (ai<a-1) ==True:
        up = scipy.linalg.block_diag(up,ut)
        down = scipy.linalg.block_diag(down,dot)
        ai+=1

    up = np.roll(up,4*nx,axis=1)
    down = np.roll(down,-4*nx,axis=1)
    
    
    a = nz
    ai = 0 
    uppertry = np.array(up)
    downtry = np.array(down)
    while (ai<a-1) ==True:
        uppertry = scipy.linalg.block_diag(uppertry,up)
        downtry = scipy.linalg.block_diag(downtry,down)
        ai+=1

    a = nx
    ai = 0 
    left = np.array(lt)
    right = np.array(rt)
    while (ai<a-1) ==True:
        left = scipy.linalg.block_diag(left,lt)
        right = scipy.linalg.block_diag(right,rt)
        ai+=1

    left = np.roll(left,-4,axis=1)
    right = np.roll(right,4,axis=1)
    
    
    a = ny*nz
    ai = 0 
    lefttry = np.array(left)
    righttry = np.array(right)
    while (ai<a-1) ==True:
        lefttry = scipy.linalg.block_diag(lefttry,left)
        righttry = scipy.linalg.block_diag(righttry,right)
        ai+=1

    
    total = diag+front+back+uppertry+downtry+lefttry+righttry
    
    beta_total = np.dot(beta_blocks,total)
    #checked and this is hermitian which is great so now can use .eigh

    eigenval,eigenvect = np.linalg.eigh(total)
    
    return beta_total,eigenval,eigenvect




