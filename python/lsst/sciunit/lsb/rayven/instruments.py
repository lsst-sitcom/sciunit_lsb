import batoid
import numpy as np

class CBP:

    def __init__(self, cbp_az, cbp_alt, rotTelPos=0):

        self.az = cbp_az
        self.alt = cbp_alt
        #self.dome_az = dome_az

        self.calculate_position(cbp_az)
        self.model = self.get_batoid_model(rotTelPos=rotTelPos)

    def calculate_position(self, angle_deg, radius=12.4, z=13.13): #radius & z in m, angle in deg

        self.x = radius * np.cos(np.deg2rad(angle_deg))
        self.y = - radius * np.sin(np.deg2rad(angle_deg)) #Azimuth is in a left handed coord system, (x,y,z) is right handed
        self.z = z

    def get_batoid_model(self, rotTelPos=0):
        cbp = batoid.Optic.fromYaml("CBP.yaml")
        # Need to rotate the fiducial model, which points upwards, to point downwards
        cbp = cbp.withLocalRotation(batoid.RotX(np.deg2rad(180)))
        # Rotate 120 degrees to make cbpAz=0 point directly at the Simonyi Survey telescope
        cbp = cbp.withLocalRotation(batoid.RotZ(np.deg2rad(120)))
        # The apply additional requested azimuth
        cbp = cbp.withLocalRotation(batoid.RotZ(np.deg2rad(self.az)))
        # Move CBP to it's position in the dome
        cbp = cbp.withGlobalShift(
            [12.4*np.cos(np.deg2rad(30)), 
             -12.4*np.sin(np.deg2rad(30)), 
             12.135 + 0.998  # Height above azimuth ring + height of CBP itself
            ]
        )
        # Now point in altitude
        cbp = cbp.withLocalRotation(batoid.RotX(np.deg2rad(90-self.alt)))
        cbp = cbp.withLocallyRotatedOptic("Cassegrain", batoid.RotZ(np.deg2rad(rotTelPos)))
        return cbp
    
class LSST:

    def __init__(self, tel_az, tel_alt, dome_az, rotTelPos=0, w_pupil=False):

        self.az = tel_az
        self.alt = tel_alt
        self.dome_az = dome_az

        self.dome_corrected_az = self.az-self.dome_az+90 #90 deg offset for the Dome

        self.model = self.get_batoid_model(rotTelPos=rotTelPos, w_pupil=w_pupil)

    def get_batoid_model(self, rotTelPos=0, w_pupil=False):
        # All in degrees
        elevBearingHeight = 5.425  # Height of elevation bearing above azimuth ring
        m1VertexHeight = 3.53 # When Zenith pointint
        
        # Fiducial model has M1 vertex at origin
        lsst = batoid.Optic.fromYaml("LSST_g.yaml")
        # Add a phantom surface to record the pupil intersection
        if w_pupil:
            print('add phantom surface')
            lsst = lsst.withInsertedOptic(
                before='M1',
                item=batoid.Interface(
                    name='pupil',
                    surface=batoid.Plane(),
                    coordSys=batoid.CoordSys(
                        origin=[0,0,0.4393899101271662]
                    )
                ),
            )    
        # Shift up so that center of azimuth ring is the origin
        lsst = lsst.withGlobalShift([0, 0, m1VertexHeight])
        # Apply camera rotation
        lsst = lsst.withLocallyRotatedOptic("LSSTCamera", batoid.RotZ(np.deg2rad(rotTelPos)))
        # Apply azimuth rotation    
        lsst = lsst.withLocalRotation(batoid.RotZ(np.deg2rad(90-self.dome_corrected_az)))  # Point in Az
        # Apply altitude rotation.  This one's a little tricky as you need to
        # specify a non-default rotation axis.
        lsst = lsst.withLocalRotation(
            batoid.RotX(np.deg2rad(90-self.alt)), 
            rotCenter=[0, 0, elevBearingHeight], 
            coordSys=batoid.globalCoordSys
        )
        return lsst
    
        
        
        
        

        
        

        
        