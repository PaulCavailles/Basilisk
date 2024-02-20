from dataclasses import dataclass

from Basilisk.fswAlgorithms import (mrpFeedback, attTrackingError,hillPoint,
                                    inertial3D, rwMotorTorque,
                                    tamComm, mtbMomentumManagementSimple, 
                                    torque2Dipole, dipoleMapping, 
                                    mtbFeedforward, rwNullSpace)

from Basilisk.utilities import SimulationBaseClass
from Basilisk.utilities import macros
from Basilisk.utilities import simSetPlanetEnvironment
from Basilisk.utilities import unitTestSupport 
from Basilisk.utilities import RigidBodyKinematics
from Basilisk.utilities import simIncludeRW
from Basilisk.utilities import (orbitalMotion,simIncludeGravBody)

from Basilisk.simulation import facetDragDynamicEffector, dragDynamicEffector
from Basilisk.simulation import magneticFieldWMM,magnetometer, MtbEffector
from Basilisk.simulation import spacecraft, simpleNav
from Basilisk.simulation import simSynch
from Basilisk.simulation import exponentialAtmosphere
from Basilisk.simulation import GravityGradientEffector
from Basilisk.simulation import starTracker
from Basilisk.simulation import reactionWheelStateEffector
from Basilisk.simulation import extForceTorque

from Basilisk.simulation import radiationPressure

from Basilisk.architecture import messaging
 
from figuresPlot.plotter import *

import numpy as np
import matplotlib.pyplot as plt

#import gravity and orbital motion


#Vizard tcp://localhost:5556
# attempt to import vizard
from Basilisk.utilities import vizSupport


#basilisk path
from Basilisk import __path__
bskPath = __path__[0]

# Plotting functions
def plot_attitude_error( timeData, dataSigmaBR):
    """Plot the attitude result."""
    plt.figure(1)
    fig = plt.gcf()
    ax = fig.gca()
    vectorData = dataSigmaBR
    sNorm = np.array([np.linalg.norm(v) for v in vectorData])
    plt.plot( timeData, sNorm,
             color=unitTestSupport.getLineColor(1, 3),
             )
    plt.xlabel('Time [min]')
    plt.ylabel(r'Attitude Error Norm $|\sigma_{B/R}|$')
    ax.set_yscale('log')

def plot_control_torque( timeData, dataLr):
    """Plot the control torque response."""
    plt.figure(2)
    for idx in range(3):
        plt.plot( timeData, dataLr[:, idx],
                 color=unitTestSupport.getLineColor(idx, 3),
                 label='$L_{r,' + str(idx) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('Control Torque $L_r$ [Nm]')

def plot_rate_error( timeData, dataOmegaBR):
    """Plot the body angular velocity tracking error."""
    plt.figure(3)
    for idx in range(3):
        plt.plot( timeData, dataOmegaBR[:, idx],
                 color=unitTestSupport.getLineColor(idx, 3),
                 label=r'$\omega_{BR,' + str(idx) + '}$')
    plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('Rate Tracking Error [rad/s] ')
    return

def plot_orientation( timeData, dataPos, dataVel, dataSigmaBN):
    """Plot the spacecraft orientation."""
    vectorPosData = dataPos
    vectorVelData = dataVel
    vectorMRPData = dataSigmaBN
    data = np.empty([len(vectorPosData), 3])
    for idx in range(0, len(vectorPosData)):
        ir = vectorPosData[idx] / np.linalg.norm(vectorPosData[idx])
        hv = np.cross(vectorPosData[idx], vectorVelData[idx])
        ih = hv / np.linalg.norm(hv)
        itheta = np.cross(ih, ir)
        dcmBN = RigidBodyKinematics.MRP2C(vectorMRPData[idx])
        data[idx] = [np.dot(ir, dcmBN[0]), np.dot(itheta, dcmBN[1]), np.dot(ih, dcmBN[2])]
    plt.figure(4)
    labelStrings = (r'$\hat\imath_r\cdot \hat b_1$'
                    , r'${\hat\imath}_{\theta}\cdot \hat b_2$'
                    , r'$\hat\imath_h\cdot \hat b_3$')
    for idx in range(3):
        plt.plot( timeData, data[:, idx],
                 color=unitTestSupport.getLineColor(idx, 3),
                 label=labelStrings[idx])
    plt.legend(loc='lower right')
    plt.xlabel('Time [min]')
    plt.ylabel('Orientation Illustration')



def run(orbit,useSphericalHarmonics,simulationTimeStep,visualisation):
    
    #------------------------------------------------------------------#
    ######################Simulation Process############################
    #------------------------------------------------------------------#
    scSim = SimulationBaseClass.SimBaseClass()
    #------------------------------------------------------------------#
    
    
    
    
    
    #------------------------------------------------------------------#
    ##############################DYNAMICS##############################
    #------------------------------------------------------------------#
    
    dynamicsTaskName = "dynamics"
    #  create the simulation dynamics process                          # 
    dynProcess = scSim.CreateNewProcess("dynamicsProcess")             #
    #create dynamics task                                              #
    dynamicsTask = scSim.CreateNewTask(dynamicsTaskName,macros.sec2nano(simulationTimeStep)) #
    #binding task to process                                           #
    dynProcess.addTask(dynamicsTask)                                   #
    #------------------------------------------------------------------#
    
   
    
    
    
    
    #------------------------------------------------------------------#
    ##############################SPACECRAFT############################
    #------------------------------------------------------------------#                          
    scObject = spacecraft.Spacecraft()                                 
    scObject.ModelTag = "Tolo-sat"
    # define the simulation inertia
    I =[0.037, 0.0, 0.0, 0.0, 0.037, 0.0, 0.0, 0.0, 0.006]
    scObject.hub.IHubPntBc_B = [[0.037, 0.0, 0.0], [0.0, 0.037, 0.0], [0.0, 0.0, 0.006]]
    scObject.hub.mHub = 2.5
    # add spacecraft object to the simulation process
    scSim.AddModelToTask(dynamicsTaskName, scObject, ModelPriority=100)
    
    #monitor spacecraft states 
    #scStates = scObject.scStateOutMsg.recorder()
    #scSim.AddModelToTask(dynamicsTaskName,scStates)
    
    #------------------------------------------------------------------#
    
    
    
    #------------------------------------------------------------------#
    ########################## NAVIGATION   ############################
    #------------------------------------------------------------------# 
    
    # add the simple Navigation sensor module.  This sets the SC attitude, rate, position
    # velocity navigation message
    sNavObject = simpleNav.SimpleNav()
    sNavObject.ModelTag = "SimpleNavigation"
    scSim.AddModelToTask(dynamicsTaskName, sNavObject,ModelPriority=5)
    
    #
    #   setup the FSW algorithm tasks
    #

    # setup hillPoint guidance module
    attGuidance = hillPoint.hillPoint()
    attGuidance.ModelTag = "hillPoint"
    attGuidance.transNavInMsg.subscribeTo(sNavObject.transOutMsg)
    # if you want to connect attGuidance.celBodyInMsg, then you need a planet ephemeris message of
    # type EphemerisMsgPayload.  In this simulation the input message is not connected to create an empty planet
    # ephemeris message which puts the earth at (0,0,0) origin with zero speed.
    CelBodyData = messaging.EphemerisMsgPayload() # make zero'd planet ephemeris message
    celBodyInMsg = messaging.EphemerisMsg().write(CelBodyData)
    attGuidance.celBodyInMsg.subscribeTo(celBodyInMsg)
    scSim.AddModelToTask(dynamicsTaskName, attGuidance)

    # setup the attitude tracking error evaluation module
    attError = attTrackingError.attTrackingError()
    attError.ModelTag = "attErrorInertial3D"
    scSim.AddModelToTask(dynamicsTaskName, attError)
    attError.sigma_R0R = [0, 1, 0]
    attError.attRefInMsg.subscribeTo(attGuidance.attRefOutMsg)
    attError.attNavInMsg.subscribeTo(sNavObject.attOutMsg)

    # setup the MRP Feedback control module
    mrpControl = mrpFeedback.mrpFeedback()
    mrpControl.ModelTag = "mrpFeedback"
    scSim.AddModelToTask(dynamicsTaskName, mrpControl)
    mrpControl.guidInMsg.subscribeTo(attError.attGuidOutMsg)
    mrpControl.K = 3.5
    mrpControl.Ki = -1.0  # make value negative to turn off integral feedback
    mrpControl.P = 30.0
    mrpControl.integralLimit = 2. / mrpControl.Ki * 0.1

    # the control torque is read in through the messaging system
    extFTObject = extForceTorque.ExtForceTorque()
    extFTObject.ModelTag = "externalDisturbance"
    # use the input flag to determine which external torque should be applied
    # Note that all variables are initialized to zero.  Thus, not setting this
    # vector would leave it's components all zero for the simulation.
    scObject.addDynamicEffector(extFTObject)
    scSim.AddModelToTask(dynamicsTaskName, extFTObject)

    # connect torque command to external torque effector
    extFTObject.cmdTorqueInMsg.subscribeTo(mrpControl.cmdTorqueOutMsg)
    
    
    #------------------------------------------------------------------#
    ########################      GRAVITY     ##########################
    #------------------------------------------------------------------# 
    
    #creation of the gravity factory
    gravFactory = simIncludeGravBody.gravBodyFactory()
    
    #create Earth
    planet = gravFactory.createEarth()
    # ensure Earth is the central gravitational body
    planet.isCentralBody = True          
    
    if useSphericalHarmonics:
        planet.useSphericalHarmParams = True
        simIncludeGravBody.loadGravFromFile(bskPath + '/supportData/LocalGravData/GGM03S-J2-only.txt',
                                            planet.spherHarm,2)
    mu = planet.mu
    
    theSun = gravFactory.createSun()
    
    # Create a sun spice message, zero it out, required by srp


    # setup spice library for Earth ephemeris and Hubble states
    timeInitString = "2015 February 10, 00:00:00.0 TDB"
    gravFactory.createSpiceInterface(bskPath + '/supportData/EphemerisData/',
                                     timeInitString,
                                     epochInMsg=True)
    gravFactory.spiceObject.zeroBase = 'Earth'
    #scNames = ["TOLOSAT"]
    #gravFactory.spiceObject.addSpacecraftNames(messaging.StringVector(scNames))
 
    
    gravFactory.spiceObject.loadSpiceKernel("hst_edited.bsp", bskPath + '/supportData/EphemerisData/')

    
    # Bind gravity model to spacecraft
    scObject.gravField.gravBodies = spacecraft.GravBodyVector(list(gravFactory.gravBodies.values()))

    # need spice to run before spacecraft module as it provides the spacecraft translational states
    
    scSim.AddModelToTask(dynamicsTaskName, gravFactory.spiceObject)
    scSim.AddModelToTask(dynamicsTaskName, scObject)
    #------------------------------------------------------------------#


   




    #------------------------------------------------------------------#
    ########################     ORBIT        ##########################
    #------------------------------------------------------------------# 
    
    # setup the orbit using classical orbit elements
    oe = orbitalMotion.ClassicElements()
    
    oe.a = orbit.a
    oe.e = orbit.e
    oe.i = orbit.i
    oe.Omega = orbit.Omega
    oe.omega = orbit.omega
    oe.f = orbit.ta
    
    
   
    
    rN, vN = orbitalMotion.elem2rv(mu, oe)
    oe = orbitalMotion.rv2elem(mu, rN, vN)
    
    
    scObject.hub.r_CN_NInit = rN  # m   - r_BN_N
    scObject.hub.v_CN_NInit = vN  # m/s - v_BN_N
    
    
    # set the simulation time
    n = np.sqrt(mu / oe.a / oe.a / oe.a) #mean anomaly
    P = 2. * np.pi / n # one orbit 

    simulationTime = macros.sec2nano(P)
    
    #------------------------------------------------------------------#
    #####################      Atmosphere        #######################
    #------------------------------------------------------------------#

    
    atmoModule = exponentialAtmosphere.ExponentialAtmosphere()
    atmoModule.ModelTag = "ExpAtmo"
    simSetPlanetEnvironment.exponentialAtmosphere(atmoModule, "earth")
    
    #atmosphere module needs spacecraft position to evaluate density 

    
    scSim.AddModelToTask(dynamicsTaskName,atmoModule, ModelPriority=3)

    
    
    #atmStates = atmoModule.envOutMsgs[0].recorder()   
    
     
    #scSim.AddModelToTask(dynamicsTaskName,atmStates)
    
    
    #------------------------------------------------------------------#
    ########################      DRAG        ##########################
    #------------------------------------------------------------------# 

  
    scDrag = facetDragDynamicEffector.FacetDragDynamicEffector()
    scDrag.ModelTag = "facetDrag"
    
    
    scAreas = [0.01, 0.01, 0.03, 0.03, 0.03, 0.03]
    scCoeff = np.array([2.2, 2.2, 2.2, 2.2, 2.2, 2.2])
    B_normals = [np.array([1, 0, 0]),
                 np.array([-1, 0, 0]),
                 np.array([0, 1, 0]),
                 np.array([0, -1, 0]),
                 np.array([0, 0, 1]),
                 np.array([0, 0, -1])]
    B_locations = [np.array([0.15,0,0]),
                   np.array([-0.15,0,0]),
                   np.array([0,0.05,0]),
                   np.array([0,-0.05,0]),
                   np.array([0.,0.,0.05]),
                   np.array([0.,0.,-0.05])]
    for ind in range(0,len(scAreas)):
       scDrag.addFacet(scAreas[ind], scCoeff[ind], B_normals[ind], B_locations[ind])
    
    # #drag module needs     
    dragEffectorTaskName = "drag"
    dynProcess.addTask(scSim.CreateNewTask(dragEffectorTaskName, macros.sec2nano(simulationTimeStep)))
    scSim.AddModelToTask(dragEffectorTaskName, scDrag)


    atmoModule.addSpacecraftToModel(scObject.scStateOutMsg)
    
    
    scObject.addDynamicEffector(scDrag)
    
    scDrag.atmoDensInMsg.subscribeTo(atmoModule.envOutMsgs[0])
    
    #scSim.AddVariableForLogging(scDrag.ModelTag+'.forceExternal_B', macros.sec2nano(1.0))
    scSim.AddVariableForLogging(scDrag.ModelTag+'.forceExternal_B', macros.sec2nano(1.0),StartIndex=0, StopIndex=2)
    
    
    #------------------------------------------------------------------#
    ########################     GRAVITY GRADIENT      #################
    #------------------------------------------------------------------# 

    ggEff = GravityGradientEffector.GravityGradientEffector()
    ggEff.ModelTag = scObject.ModelTag
    ggEff.addPlanetName(planet.planetName)
    scObject.addDynamicEffector(ggEff)
    scSim.AddModelToTask(dynamicsTaskName, ggEff)
  
    #------------------------------------------------------------------#
    ########################   Magnetic Field  #########################
    #------------------------------------------------------------------#
    
    # create magnetic field module
    magModule = magneticFieldWMM.MagneticFieldWMM()
    magModule.ModelTag = "WMM"
    magModule.dataPath = bskPath + '/supportData/MagneticField/'
    epochMsg = unitTestSupport.timeStringToGregorianUTCMsg('2019 June 27, 10:23:0.0 (UTC)')
    magModule.epochInMsg.subscribeTo(epochMsg)
    magModule.addSpacecraftToModel(scObject.scStateOutMsg)  # this command can be repeated if multiple
    scSim.AddModelToTask(dynamicsTaskName, magModule)
    
    
    # add magnetic torque bar effector
    # mtbEff = MtbEffector.MtbEffector()
    # mtbEff.ModelTag = "MtbEff"
    # scObject.addDynamicEffector(mtbEff)
    # scSim.AddModelToTask(dynamicsTaskName, mtbEff)    
    
    # create the minimal TAM module
    TAM = magnetometer.Magnetometer()
    TAM.ModelTag = "TAM_sensor"
    # specify the optional TAM variables
    TAM.scaleFactor = 1.0
    TAM.senNoiseStd = [0.000001,  0.000001, 0.000001]
    scSim.AddModelToTask(dynamicsTaskName, TAM)
    
    # setup tamComm module
    tamCommConfig = tamComm.tamConfigData()
    tamCommConfig.dcm_BS = [1., 0., 0., 0., 1., 0., 0., 0., 1.]
    tamCommWrap = scSim.setModelDataWrap(tamCommConfig)
    tamCommWrap.ModelTag = "tamComm"
    scSim.AddModelToTask(dynamicsTaskName, tamCommWrap, tamCommConfig)

    
    #------------------------------------------------------------------#
    ########################    SRP      ###############################
    #------------------------------------------------------------------#
    
    # Create an SRP model
    srp = radiationPressure.RadiationPressure()  # default model is the SRP_CANNONBALL_MODEL
    srp.area = 0.1*0.01 + 0.3*0.1*4  # m^3
    srp.coefficientReflection = 0.6
    srp.ModelTag="srp_module"
    scSim.AddModelToTask(dynamicsTaskName, srp, ModelPriority=90)
    scObject.addDynamicEffector(srp)
    
    scSim.AddVariableForLogging(srp.ModelTag+".forceExternal_N",macros.sec2nano(1.),StartIndex=0, StopIndex=2  )
        
    #------------------------------------------------------------------#
    ########################    R W      ###############################
    #------------------------------------------------------------------# 
    
    # useJitterSimple = False
    # useRWVoltageIO = False
    
    # varRWModel = messaging.BalancedWheels
    # if useJitterSimple:
    #     varRWModel = messaging.JitterSimple
    
    # rwFactory = simIncludeRW.rwFactory()
    # RW = rwFactory.create('Honeywell_HR16', [1,0,0], maxMomentum = 50.0, Omega = 100.0, RWModel = varRWModel )    
    # numRW = rwFactory.getNumOfDevices()
    # rwStateEffector = reactionWheelStateEffector.ReactionWheelStateEffector()
    # rwStateEffector.ModelTag = "RWcluster"
    # rwFactory.addToSpacecraft(scObject.ModelTag, rwStateEffector, scObject)

    # # add RW object array to the simulation process.  This is required for the UpdateState() method
    # # to be called which logs the RW states
    # scSim.AddModelToTask(dynamicsTaskName, rwStateEffector, None, 20)
    
    #------------------------------------------------------------------#
    ###########################VIZARD BIND #############################
    #------------------------------------------------------------------#
    if vizSupport.vizFound:
        viz = vizSupport.enableUnityVisualization(scSim,"dynamics",
                                                  scObject, liveStream=visualisation
                                                  # , saveFile=fileName
                                                  )
        viz.settings.spacecraftSizeMultiplier = 2.5
        viz.settings.mainCameraTarget='earth'
      
    #------------------------------------------------------------------#
    
    
    
    
    
    #------------------------------------------------------------------#
    ###########################SIMULATION START#########################
    #------------------------------------------------------------------#
    
    # Setup data logging before the simulation is initialized
    dataRec = scObject.scStateOutMsg.recorder(macros.sec2nano(1.))
    ggLog = ggEff.gravityGradientOutMsg.recorder(macros.sec2nano(1.))
    tamLog = TAM.tamDataOutMsg.recorder(macros.sec2nano(1.))
    tamCommLog = tamCommConfig.tamOutMsg.recorder(macros.sec2nano(1.))
    magLog = magModule.envOutMsgs[0].recorder(macros.sec2nano(1.))
    spiceDataLog = gravFactory.spiceObject.planetStateOutMsgs[1].recorder(macros.sec2nano(1.))
    attErrLog = attError.attGuidOutMsg.recorder(macros.sec2nano(1.))
    mrpLog = mrpControl.cmdTorqueOutMsg.recorder(macros.sec2nano(1.))
    snAttLog = sNavObject.attOutMsg.recorder(macros.sec2nano(1.))
    snTransLog = sNavObject.transOutMsg.recorder(macros.sec2nano(1.))


    #dataAtmoLog = atmoModule.envOutMsgs[0].recorder(1)
    scSim.AddModelToTask(dynamicsTaskName , dataRec)
    scSim.AddModelToTask(dynamicsTaskName , tamLog)
    scSim.AddModelToTask(dynamicsTaskName , tamCommLog)
    scSim.AddModelToTask(dynamicsTaskName, magLog)
    scSim.AddModelToTask(dynamicsTaskName,ggLog)
    scSim.AddModelToTask(dynamicsTaskName,spiceDataLog)
    scSim.AddModelToTask(dynamicsTaskName, mrpLog)
    scSim.AddModelToTask(dynamicsTaskName, attErrLog)
    scSim.AddModelToTask(dynamicsTaskName, snAttLog)
    scSim.AddModelToTask(dynamicsTaskName, snTransLog)

    # create the FSW vehicle configuration message
    vehicleConfigOut = messaging.VehicleConfigMsgPayload()
    vehicleConfigOut.ISCPntB_B = I  # use the same inertia in the FSW algorithm as in the simulation
    configDataMsg = messaging.VehicleConfigMsg().write(vehicleConfigOut)
    mrpControl.vehConfigInMsg.subscribeTo(configDataMsg)


    
    
    #scSim.AddModelToTask(dynamicsTaskName, dataAtmoLog)
    
    
    #Messages connections
    
    
    sNavObject.scStateInMsg.subscribeTo(scObject.scStateOutMsg)
    TAM.stateInMsg.subscribeTo(scObject.scStateOutMsg)
    TAM.magInMsg.subscribeTo(magModule.envOutMsgs[0])
    tamCommConfig.tamInMsg.subscribeTo(TAM.tamDataOutMsg)
    srp.sunEphmInMsg.subscribeTo(gravFactory.spiceObject.planetStateOutMsgs[1])
    
    scSim.ShowExecutionOrder()
    scSim.InitializeSimulation()
    scSim.ConfigureStopTime(simulationTime)
    
    scSim.ExecuteSimulation()
    

    timeData = dataRec.times() * macros.NANO2MIN
    dragBodyData = scSim.GetLogVariableData(scDrag.ModelTag+".forceExternal_B")
    srpForceData = scSim.GetLogVariableData(srp.ModelTag+".forceExternal_N")
    dataLr = mrpLog.torqueRequestBody
    dataSigmaBR = attErrLog.sigma_BR
    dataOmegaBR = attErrLog.omega_BR_B
    dataPos = snTransLog.r_BN_N
    dataVel = snTransLog.v_BN_N
    dataSigmaBN = snAttLog.sigma_BN
    
    print("sunName = ",theSun.planetName)
    print(ggLog.gravityGradientTorque_B.shape)
    num = 1
    plotSCpos(timeData, dataRec,num)
    num +=1
    plotEarthMagneticField(timeData, magLog,num)
    num +=1
    plotSensorMagneticField(timeData, tamLog, num)
    num +=1 
    plotGGbodyFrame(timeData, ggLog,num)
    num+=1
    plotDragEffector(timeData, dragBodyData,num)
    num+=1
    plotSpicePlanetData(timeData,spiceDataLog,num)
    num+=1
    plotSRPEffector(timeData,srpForceData,num)

    plot_attitude_error( timeData, dataSigmaBR)
    num+=1

    plot_control_torque( timeData, dataLr)
    
    plot_rate_error( timeData, dataOmegaBR)
    num+=1
    plot_orientation( timeData, dataPos, dataVel, dataSigmaBN) 
    num+=1
    
    scSim.ShowExecutionFigure(True)
    
    

############################## Orbit Structure ##################
class Orbit:    
    def __init__(self,a , e, i, Omega, omega, ta):
        self.a = a
        self.e = e
        self.i = i
        self.Omega = Omega
        self.omega = omega
        self.ta = ta 







if __name__ == "__main__":
    
    orbit = Orbit(a = (6371 + 400) * 1000,
                  e = 0.002,
                  i= 97.4 * macros.D2R,
                  Omega= 10.70 * macros.D2R,
                  omega= 90 * macros.D2R,
                  ta = 0.1 * macros.D2R
                )
    
    run(
        orbit,
        True,# useSphericalHarmonics
        1.0, # time step (s)
        False # liveStream
    )
