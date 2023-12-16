import numpy as np
import matplotlib.pyplot as plt

def plotSCpos(timeData, dataRec,num):
    plt.figure(num)
    dim = 3
    for idx in range(dim):
        plt.plot(timeData,dataRec.r_BN_N[:,idx]/1000.0)
    
    plt.xlabel("Time [mn]")
    plt.ylabel("Position [Km]")
    plt.title("Satellite position w.r.t inertial")
    plt.grid()
    
    
def plotEarthMagneticField(timeData, dataLog,num):
    plt.figure(num)
    dim = 3
    plt.plot(timeData,dataLog.magField_N[:,0]* 1e9,
             timeData,dataLog.magField_N[:,1]* 1e9,
             timeData,dataLog.magField_N[:,2]* 1e9)
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.xlabel("Time [min]")
    plt.ylabel("Magnetic Field [nT]")
    plt.title("Earth magnetique field seen by the satellite")
    plt.grid()
    
def plotSensorMagneticField(timeData, dataLog, num):
    plt.figure(num)
    plt.plot(timeData,dataLog.tam_S[:,0]* 1e9,
             timeData,dataLog.tam_S[:,1]* 1e9,
             timeData,dataLog.tam_S[:,2]* 1e9)
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.xlabel("Time [min]")
    plt.ylabel("Magnetic Field [nT]")
    plt.title("Magnetic field seen by the sensor")
    plt.grid()
    
def plotGGbodyFrame(timeData, dataLog,num):
     #gravityGradientTorque_B
    plt.figure(num)
    plt.plot(timeData,dataLog.gravityGradientTorque_B[:,0],
             timeData,dataLog.gravityGradientTorque_B[:,1],
             timeData,dataLog.gravityGradientTorque_B[:,2])
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.xlabel("Time [min]")
    plt.ylabel("Gravity Gradient [Nm]")
    plt.title("Gravity Gradient in body frame")
    plt.grid()
    
def plotDragEffector(timeData, dataLog,num):
     #gravityGradientTorque_B
    plt.figure(num)
    plt.plot(timeData,dataLog[:,1],
             timeData,dataLog[:,2],
             timeData,dataLog[:,3])
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.legend(["x", "y", "z" ])
    plt.xlabel("Time [min]")
    plt.ylabel("Drag on each face [N]")
    plt.title("Drag in body frame")
    plt.grid()



def plotSRPEffector(timeData, dataLog,num):
    plt.figure(num)
    plt.plot(timeData,dataLog[:,1],
             timeData,dataLog[:,2],
             timeData,dataLog[:,3])
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.legend(["x", "y", "z" ])
    plt.xlabel("Time [min]")
    plt.ylabel("SRP in each direction [N]")
    plt.title("SRP force in inertial frame")
    plt.grid()
    
def plotSpicePlanetData(timeData, dataLog,num):
    plt.figure(num)
    dim = 3
    plt.plot(timeData,dataLog.PositionVector[:,0],
             timeData,dataLog.PositionVector[:,1],
             timeData,dataLog.PositionVector[:,2],)
    plt.autoscale(enable=True, axis='both', tight=None)
    plt.xlabel("Time [min]")
    plt.ylabel("position [meters]")
    plt.title("Earth ephemeris")
    plt.legend(['x','y','z'])
    plt.grid()