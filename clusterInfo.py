#This python code is for ORELSE Cluster information setup
import numpy as np

def clusterInfo(name)

    """
    Note:
    input:
        name    ---- cluster
    output:
        zmax    ---- MAX redshift for field
        zmin    ---- MIN redshift for field
        label   ---- label of cluster's name, could be used in position plot
        circles ---- RA & Dec for cluster/groups center
        centerZ ---- redshift for each cluster
        vel     ---- velocity dispersion for each cluster in 0.5Mpc [km/s]
        velErr  ---- uncertainty of velocity disoersion
        deltaZ  ---- uncertainty of redshift for the cluster/groups
    """
    ###Cl0023
    if (name == 'Cl0023'):
        zmax = 0.82
        zmin = 0.87
        label = ['A', 'B1','B2', 'C', 'M']
        circles = ((6.02560,4.3590),(5.97570,4.3884),(5.96970,4.3820),(5.92470,4.3807),(5.96740,4.3199))
        centerZ = (0.8396, 0.8290, 0.8453, 0.8466, 0.8472)
        vel = (507.0, 106.2, 231.3, 543.8, 487.3)
        velErr = (125.5, 51.4, 53.8, 58.7, 84.8)

    ###Cl1604
    if (Name == 'Cl1604'):
        zmax = 0.84
        zmin = 0.96
        label = ['A', 'B', 'C', 'D', 'F','G','H','I']
        #C G H I are groups
        centerZ = (0.8984, 0.8648, 0.9344, 0.9277, 0.9331, 0.9019, 0.8528, 0.9024)
        circles = ((241.09311,43.0821),(241.10796,43.2397),(241.03142,43.2679),(241.14094,43.3539),
                   (241.20104,43.3684),(240.92745,43.4030),(240.89890,43.3669),(240.79746,43.3915))
        vel = (576.9, 792.9, 1079.6, 675.6, 619.9, 398.5, 283.2, 163.0)
        velErr = (120.8, 80.1, 289.6, 180.1, 135.0, 85.4, 71.7, 65.1)

    ###Cl1324
    if (Name == 'Cl1324'):
        zmax = 0.65
        zmin = 0.79
        label = ['A', 'B', 'C', 'I']
        #C is a group
        centerZ = (0.7556, 0.6979, 0.7382, 0.6957)
        circles = ((201.20097, 30.1920),(201.08830, 30.2156),(201.00770, 30.4164),(201.20096,30.9680))
        vel = (1019.2, 897.0, 143.8, 646.1)
        velErr = (142.0, 394.7, 40.7, 113.0)


    ###NEP5281
    if (name == 'NEP5281')
        zmax = 0.80
        zmin = 0.84
        label = ['A', 'Null']
        circles = ((275.38451, 68.465768),(0.,0.))
        centerZ = (0.8168, 0)
        vel = (1146.1, 0.)
        velErr = (124.8, 0.)


    ###NEP200
    if (name == 'NEP200')
        zmax = 0.68
        zmin = 0.71
        label = ['A', 'Null']
        circles = ((269.33196, 66.525991),(0.,0.))
        centerZ = (0.6931, 0)
        vel = (541.1, 0.)
        velErr = (139.2, 0.)

    c = 3 * 10 **5 #[km/s]
    deltaZ = 3 * np.divide(vel, c)

    return zmax, zmin, label, circles, centerZ, vel, velErr, deltaZ