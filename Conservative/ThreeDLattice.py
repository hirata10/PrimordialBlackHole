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

#can get dis_rel_xyz without wilson if rx=ry=rz=0 i think  
def disp_rel_xyzW(nx,ny,nz,Ex,Ey,Ez,m,rx,ry,rz,x_flag=True,y_flag=True,z_flag=True):
    
    dx = Ex
    dy = Ey
    dz = Ez
    kx= np.zeros(nx)
    ky= np.zeros(ny)
    kz= np.zeros(nz)
    
    if x_flag ==False and y_flag==False and z_flag==False:
        print('not valid. Must fourier transform in at least one direction')
        return None,None,None,None
    if x_flag ==True and y_flag==True and z_flag==True:
        
        kx = np.linspace(-np.pi/dx,np.pi/dx, num=nx, axis=0)

        ky = np.linspace(-np.pi/dy,np.pi/dy, num=ny, axis=0)
        kz = np.linspace(-np.pi/dz,np.pi/dz, num=nz, axis=0)
        


        diag = -m*np.identity(4) + (rx+ry+rz)*np.identity(4)
        print(diag)
        front = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)

        testvals = []
        kmags = []
        i = 0 
        
    if x_flag ==False and y_flag==True and z_flag==True:
        ky = np.linspace(-np.pi/dy,np.pi/dy, num=ny, axis=0)
        kz = np.linspace(-np.pi/dz,np.pi/dz, num=nz, axis=0)


        diag = -m*np.identity(4*nx) + (rx+ry+rz)*np.identity(4*nx)
        print(diag)
        front1 = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up1 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
        
        front = np.array(front1)
        up = np.array(up1)
        right = np.array(right1)
        for i in range(1,nx,1):
            front = scipy.linalg.block_diag(front,front1)
            up =  scipy.linalg.block_diag(up,up1)
            right =  scipy.linalg.block_diag(right,right1)
            
        
        testvals = []
        kmags = []
        i = 0 
        right = np.roll(right,4,axis=1)
    if x_flag ==True and y_flag==False and z_flag==True:
        kx = np.linspace(-np.pi/dx,np.pi/dx, num=nx, axis=0)
        kz = np.linspace(-np.pi/dz,np.pi/dz, num=nz, axis=0)


        diag = -m*np.identity(4*ny) + (rx+ry+rz)*np.identity(4*ny)
        print(diag)
        front1 = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up1 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
        
        front = np.array(front1)
        up = np.array(up1)
        right = np.array(right1)
        for i in range(1,ny,1):
            front = scipy.linalg.block_diag(front,front1)
            up =  scipy.linalg.block_diag(up,up1)
            right =  scipy.linalg.block_diag(right,right1)
        testvals = []
        kmags = []
        i = 0 
          
        up = np.roll(up,4,axis=1)
    if x_flag ==True and y_flag==True and z_flag==False:
        kx = np.linspace(-np.pi/dx,np.pi/dx, num=nx, axis=0)
        ky = np.linspace(-np.pi/dy,np.pi/dy, num=ny, axis=0)


        diag = -m*np.identity(4*nz) + (rx+ry+rz)*np.identity(4*nz)
        print(diag)
        front1 = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up1 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
        
        front = np.array(front1)
        up = np.array(up1)
        right = np.array(right1)
        for i in range(1,nz,1):
            front = scipy.linalg.block_diag(front,front1)
            up =  scipy.linalg.block_diag(up,up1)
            right =  scipy.linalg.block_diag(right,right1)
        testvals = []
        kmags = []
        i = 0 
          
        front = np.roll(front,4,axis=1)
        
    if x_flag ==True and y_flag==False and z_flag==False:
        kx = np.linspace(-np.pi/dx,np.pi/dx, num=nx, axis=0)
        


        diag = -m*np.identity(4*nz*ny) + (rx+ry+rz)*np.identity(4*nz*ny)
        print(diag)
        front1 = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up1 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
        
        front = np.array(front1)
        up = np.array(up1)
        right = np.array(right1)
        for i in range(1,nz*ny,1):
            front = scipy.linalg.block_diag(front,front1)
            up =  scipy.linalg.block_diag(up,up1)
            right =  scipy.linalg.block_diag(right,right1)
        testvals = []
        kmags = []
        i = 0 
          
        front = np.roll(front,4,axis=1)
        up = np.roll(up,4*nz,axis=1)   
    
    if x_flag ==False and y_flag==True and z_flag==False:
        ky = np.linspace(-np.pi/dy,np.pi/dy, num=ny, axis=0)
        


        diag = -m*np.identity(4*nz*nx) + (rx+ry+rz)*np.identity(4*nz*nx)
        print(diag)
        front1 = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up1 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
        
        front = np.array(front1)
        up = np.array(up1)
        right = np.array(right1)
        for i in range(1,nz*nx,1):
            front = scipy.linalg.block_diag(front,front1)
            up =  scipy.linalg.block_diag(up,up1)
            right =  scipy.linalg.block_diag(right,right1)
        testvals = []
        kmags = []
        i = 0 
          
        front = np.roll(front,4,axis=1)
        right = np.roll(right,4*nz,axis=1)   
    if x_flag ==False and y_flag==False and z_flag==True:
        kz = np.linspace(-np.pi/dz,np.pi/dz, num=nz, axis=0)
        


        diag = -m*np.identity(4*nz*nx) + (rx+ry+rz)*np.identity(4*ny*nx)
        print(diag)
        front1 = np.array( [[0,0,-1.,0],[0,0,0,1.],[-1.,0,0,0],[0,1.,0,0]])/(2*Ez) - rz*.5*np.identity(4)
        up1 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])/(2*Ey) - ry*.5*np.identity(4)
        right1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])/(2*Ex) - rx*.5*np.identity(4)
        
        front = np.array(front1)
        up = np.array(up1)
        right = np.array(right1)
        for i in range(1,ny*nx,1):
            front = scipy.linalg.block_diag(front,front1)
            up =  scipy.linalg.block_diag(up,up1)
            right =  scipy.linalg.block_diag(right,right1)
        testvals = []
        kmags = []
        i = 0 
          
        up = np.roll(front,4,axis=1)
        right = np.roll(right,4*ny,axis=1)   
    
    
    for x in range(nx):
        for y in range(ny):
            for z in range(nz):
                matr = diag + front*np.exp(complex(0,1)*Ez*kz[z]) + np.conjugate(front.T)*np.exp(complex(0,1)*-Ez*kz[z]) + up*np.exp(complex(0,1)*Ey*ky[y]) + np.conjugate(up.T)*np.exp(complex(0,1)*-Ey*ky[y])+ right*np.exp(complex(0,1)*Ex*kx[x]) + np.conjugate(right.T)*np.exp(complex(0,1)*-Ex*kx[x])
                test1val, test1vect = np.linalg.eigh(matr)
                #print(test1val)
                #kmags.append(np.sqrt(kx[x]**2 + ky[y]**2 + kz[z]**2))
                testvals.append(test1val)
                #print(i,x,y,z)
                i+=1


    return kx,ky,kz,testvals 


def disp_rel_rthetaphi(nr,ntheta,nphi,Er,Etheta,Ephi,m,rr,rtheta,rphi,rs,theta_points,phi_flag=True):
    
    dr = Er
    dtheta = Etheta
    dphi = Ephi
    

    if phi_flag==True:
        kphi = np.linspace(-np.pi/dphi,np.pi/dphi, num=nphi, axis=0)
    if phi_flag==False:
        kphi= np.zeros(nphi)
        
        
        
    
    rs_add = np.array([[0,0,0,-1],[0,0,-1,0],[0,-1,0,0],[-1,0,0,0]])/2/Er
    #need to rework the fourier transform by hand first to see whats up 
    
    
    testvals = []
    kmags = []
    i = 0 
    rs_add_1 = np.array([[0,0,0,-1],[0,0,-1,0],[0,-1,0,0],[-1,0,0,0]])/2/Er
    #rs_add = scipy.linalg.block_diag(rs_add)
    rs_add = np.array(rs_add_1)
    theta_add_1 = np.array([[0,0,0,complex(0,1)],[0,0,complex(0,-1),0],[0,complex(0,1),0,0],[complex(0,-1),0,0,0]])*np.sqrt(1-2*M/rs[0])/rs[0]/2/Etheta
    theta_add = np.array(theta_add_1)
    diag_1 = -m*np.identity(4)*np.sqrt(1-2*M/rs[0])  
    diag = np.array(diag_1) 
    phi_add_1 = phi_add = np.array([[0,0,1.,0],[0,0,0,-1.],[1.,0,0,0],[0,-1.,0,0]])*np.sqrt(1-2*M/rs[0])/rs[0]/np.sin(theta_points[0])/2/Ephi
    phi_add = np.array(phi_add_1)
    
    for i in range(0,nr,1):
        diag_1 = -m*np.identity(4)*np.sqrt(1-2*M/rs[i])  
        theta_add_1 = np.array([[0,0,0,complex(0,1)],[0,0,complex(0,-1),0],[0,complex(0,1),0,0],[complex(0,-1),0,0,0]])*np.sqrt(1-2*M/rs[i])/rs[i]/2/Etheta
        for theta_index in range(0,ntheta,1):
            if i>0 or theta_index>0:
                rs_add = scipy.linalg.block_diag(rs_add,rs_add_1)
                diag = scipy.linalg.block_diag(diag,diag_1)
                theta_add = scipy.linalg.block_diag(theta_add,theta_add_1)
                phi_add_1 = np.array([[0,0,1.,0],[0,0,0,-1.],[1.,0,0,0],[0,-1.,0,0]])*np.sqrt(1-2*M/rs[i])/rs[i]/np.sin(theta_points[theta_index])/2/Ephi
                phi_add = scipy.linalg.block_diag(phi_add,phi_add_1)


    print(rs_add.shape)
    print(phi_add.shape)
    print(diag.shape)
    
    for phi_index in range(0,nphi,1):
            #fix sizing here 

            matr = diag + np.roll(rs_add,4*ntheta,axis=1) + np.conjugate(np.roll(rs_add,-4*ntheta,axis=1).T)+ np.roll(theta_add,4,axis=1)+ np.conjugate(np.roll(theta_add,-4,axis=1).T)+phi_add*np.exp(complex(0,1)*Ephi*kphi[phi_index]) + np.conjugate(phi_add.T)*np.exp(complex(0,1)*-Ephi*kphi[phi_index])
            test1val, test1vect = np.linalg.eigh(matr)
            #print(test1val)
            print(test1vect.shape)
            testvals.append(test1val)
            #print(i,x,y,z)
            i+=1

                    
    return kphi,testvals 


def disp_rel_rthetaphiW(nr,ntheta,nphi,Er,Etheta,Ephi,m,rr,rtheta,rphi,rs,theta_points,r_flag=False,theta_flag=False,phi_flag=True):

  #10/17 Emily: If i am not mistaked, I think this part needed some work. I will take a look and see whether I am misremembering or not. 
    
    dr = Er
    dtheta = Etheta
    dphi = Ephi
    
    if r_flag==True:
        kr = np.linspace(-np.pi/dr,np.pi/dr, num=nr, axis=0)
    if r_flag==False:
        kr = np.zeros(nr)
                  
    if theta_flag==True:
        ktheta = np.linspace(-np.pi/dtheta,np.pi/dtheta, num=ntheta, axis=0)
    if theta_flag==False:
        ktheta= np.zeros(ntheta)
    
    if phi_flag==True:
        kphi = np.linspace(-np.pi/dphi,np.pi/dphi, num=nphi, axis=0)
    if phi_flag==False:
        kphi= np.zeros(nphi)
        
        
        
    
    rs_add = np.array([[0,0,0,-1],[0,0,-1,0],[0,-1,0,0],[-1,0,0,0]])/2/Er
    #need to rework the fourier transform by hand first to see whats up 
    
    
    gamma1 = np.array([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]])
    gamma2 = np.array([[0,0,0,complex(0,-1)],[0,0,complex(0,1),0],[0,complex(0,-1),0,0],[complex(0,1),0,0,0]])
    gamma3 = np.array([[0,0,1.,0],[0,0,0,-1.],[1.,0,0,0],[0,-1.,0,0]])

    #Gamma_theta = 1/2 np.sqrt((1-2M/r))**-1 *gamma1*gamma2
    #Gamma_phi = 1/2 (sin(theta)*np.sqrt((1-2M/r))**-1 *gamma1*gamma3 + cos(theta)*gamm2*gamma3)


    testvals = []
    kmags = []
    i = 0 
    rs_add = np.array([[0,0,0,-1],[0,0,-1,0],[0,-1,0,0],[-1,0,0,0]])/2/Er
    for r in range(0,nr,1):
        diag = -m*np.identity(4)*np.sqrt(1-2*M/rs[r])   
        theta_add = np.array([[0,0,0,complex(0,1)],[0,0,complex(0,-1),0],[0,complex(0,1),0,0],[complex(0,-1),0,0,0]])*np.sqrt(1-2*M/r_points[r])/r_points[r]/2/Etheta
        g_rr = 1/(1-2*M/rs[r])
        g_thetatheta = rs[r]**2
        Gamma_theta = (np.sqrt((1-2*M/rs[r]))**-1) *gamma1*gamma2/2
        for theta in range(0,ntheta,1):
            g_phiphi = (rs[r]*np.sin(theta_points[theta]))**2
        
            Gamma_phi =(np.sin(theta_points[theta])*(np.sqrt(1-2*M/rs[r])**-1) *gamma1*gamma3 + np.cos(theta_points[theta])*gamma2*gamma3)/2
        
            phi_add = np.array([[0,0,1.,0],[0,0,0,-1.],[1.,0,0,0],[0,-1.,0,0]])*np.sqrt(1-2*M/rs[r])/rs[r]/np.sin(theta_points[theta])/2/Ephi
            for phi in range(0,nphi,1):
                sqrt_neg_g = np.sqrt(2*(1-2*M/rs[r])+rs[r]**2+(rs[r]*np.sin(theta_points[theta]))**2)
                
                diag_W = sqrt_neg_g*(g_rr*Etheta*Ephi + g_thetatheta*Er*Ephi+g_phiphi*Er*Etheta)
                
                Rupup =  np.sqrt(2*(1-2*M/rs[r-1])+rs[r-1]**2+(rs[r-1]*np.sin(theta_points[theta]))**2)*1/(1-2*M/rs[r-1])*Etheta*Ephi
                Thetaupup = sqrt_neg_g*g_thetatheta*(Er*Ephi)*(1+Etheta*Gamma_theta)
                Phiupup = sqrt_neg_g*g_phiphi*(Er*Etheta)*(1+Ephi*Gamma_phi)
                
                Phiplus = sqrt_neg_g*g_phiphi*(Er*Etheta)*(1-Ephi*Gamma_phi)*(-1)*np.exp(complex(0,1)*kphi[phi]*Ephi) #eval at i
                Phiminus = sqrt_neg_g*g_phiphi*(Er*Etheta)*(1-Ephi*Gamma_phi)*(-1)*np.exp(complex(0,-1)*kphi[phi]*Ephi) #eval at i - phi
                
                Thetaplus = sqrt_neg_g*g_thetatheta*(Er*Ephi)*(1-Etheta*Gamma_theta)*(-1)*np.exp(complex(0,1)*ktheta[theta]*Etheta)#eval at i
                Thetaminus = np.sqrt(2*(1-2*M/rs[r])+rs[r]**2+(rs[r]*np.sin(theta_points[theta-1]))**2)*g_thetatheta*(Er*Ephi)*(1-Etheta*Gamma_theta)*(-1)*np.exp(complex(0,-1)*ktheta[theta]*Etheta)#eval at i - theta
                
                Rplus = sqrt_neg_g*(g_rr*Etheta*Ephi*-1)*np.exp(complex(0,1)*kr[r]*Er) #eval at i
                Rminus =  np.sqrt(2*(1-2*M/rs[r-1])+rs[r-1]**2+(rs[r-1]*np.sin(theta_points[theta]))**2)*1/(1-2*M/rs[r-1])*Etheta*Ephi*-1*np.exp(complex(0,-1)*kr[r]*Er) #eval at i - r
                
                matr = diag + rs_add*np.exp(complex(0,1)*Er*kr[r]) + np.conjugate(rs_add.T)*np.exp(complex(0,1)*-Er*kr[r]) + theta_add*np.exp(complex(0,1)*Etheta*ktheta[theta]) + np.conjugate(theta_add.T)*np.exp(complex(0,1)*-Etheta*ktheta[theta])+ phi_add*np.exp(complex(0,1)*Ephi*kphi[phi]) + np.conjugate(phi_add.T)*np.exp(complex(0,1)*-Ephi*kphi[phi])
                matr = matr + diag_W + Rupup+Thetaupup+Phiupup+ Phiplus + Phiminus + Thetaplus+Thetaminus + Rplus+Rminus
                test1val, test1vect = np.linalg.eigh(matr)
                #print(test1val)
                kmags.append(np.sqrt(kr[r]**2 + ktheta[theta]**2 + kphi[phi]**2))
                testvals.append(test1val)
                #print(i,x,y,z)
                i+=1

                    
    return kr,ktheta,kphi,testvals 
