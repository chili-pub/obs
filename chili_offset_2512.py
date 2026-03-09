from astropy import units as u
from astropy.coordinates import SkyCoord
import pdb


x_guider1 = 318 * u.arcsec  
y_guider1 = 240 * u.arcsec
pa_guider1 = 0         # IFS导星视场朝北

x_guider2 = 648 * u.arcsec
y_guider2 = 442.2 * u.arcsec
pa_guider2 = 132         # 耐焦导星视场的上方向和北的夹角 

x_ifs = 65 * u.arcsec       
y_ifs = 71 * u.arcsec
pa_ifs = pa_guider1

distance_guider1_ifs = (651.563 - 30) * u.arcsec


def gold2ifs(ra_guider, dec_guider, rot_angle):

    c_guider = SkyCoord(ra_guider, dec_guider, unit=(u.hourangle, u.deg))

    position_angle = (rot_angle + pa_guider1) * u.deg + 180 * u.deg
    separation = distance_guider1_ifs

    c_ifs = c_guider.directional_offset_by(position_angle, separation)  

    ra_ifs = ra_string_format(c_ifs.ra.to_string(u.hour))
    dec_ifs = dec_string_format(c_ifs.dec.to_string(u.degree))

    return ra_ifs, dec_ifs


def gnew2gold(ra_guider2, dec_guider2, rot_angle):
    
    c_guider2 = SkyCoord(ra_guider2, dec_guider2, unit=(u.hourangle, u.deg))

    # position_angle = 271.07611342 * u.deg 
    # separation = 478.669068 * u.arcsec
    position_angle = (271.07611342 + rot_angle) * u.deg 
    separation = 494.574 * u.arcsec

    c_guider1 = c_guider2.directional_offset_by(position_angle, separation)  
    
    ra_guider1 = ra_string_format(c_guider1.ra.to_string(u.hour))
    dec_guider1 = dec_string_format(c_guider1.dec.to_string(u.degree))

    return ra_guider1, dec_guider1



def offsets(ra1, dec1, ra2, dec2):
    
    pos1 = SkyCoord(ra1, dec1, frame='icrs', unit=(u.hourangle, u.deg))
    pos2 = SkyCoord(ra2, dec2, frame='icrs', unit=(u.hourangle, u.deg))
    dra, ddec = pos1.spherical_offsets_to(pos2)
    print('delta_ra: ', dra.to(u.arcsec))
    print('delta_dec: ', ddec.to(u.arcsec))


def ra_string_format(ra):
    
    tmp = ra.split('h')
    ra1 = tmp[0]
    tmptmp = tmp[1].split('m')
    ra2 = tmptmp[0]
    tmptmptmp = tmptmp[1].split('s')
    ra3 = tmptmptmp[0]
    
    return ra1+':'+ra2+':'+ra3

    
def dec_string_format(dec):
    
    tmp = dec.split('d')
    dec1 = tmp[0]
    tmptmp = tmp[1].split('m')
    dec2 = tmptmp[0]
    tmptmptmp = tmptmp[1].split('s')
    dec3 = tmptmptmp[0]
    
    return dec1+':'+dec2+':'+dec3



def offset4(ra_gnew, dec_gnew, ra_target, dec_target, rot_angle):
    
    ra_guider1, dec_guider1 = gnew2gold(ra_gnew, dec_gnew, rot_angle)
    
    ra_ifs, dec_ifs = gold2ifs(ra_guider1, dec_guider1, rot_angle)
    
    offsets(ra_ifs, dec_ifs, ra_target, dec_target)
    # print(ra_ifs, dec_ifs, ra_target, dec_target)

def gnew2others(ra_gnew, dec_gnew, rot_angle, outfile):

    ra_guider1, dec_guider1 = gnew2gold(ra_gnew, dec_gnew, rot_angle)
    
    ra_ifs, dec_ifs = gold2ifs(ra_guider1, dec_guider1, rot_angle)

    write_reg(outfile, ra_gnew, dec_gnew, 
              ra_guider1, dec_guider1, 
              ra_ifs, dec_ifs, rot_angle)

def write_reg(outfile, 
              ra_guider2, dec_guider2, 
              ra_guider1, dec_guider1, 
              ra_ifs, dec_ifs, rot_angle):
    
    with open(outfile, "w") as f:
        f.write('# Region file format: DS9 version 4.1 \n')
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        f.write('fk5 \n')
        
        str_guider2 = str_box(ra_guider2, dec_guider2, 
                              x_guider2.value, y_guider2.value, 
                              angle=pa_guider2+rot_angle)
        f.write(str_guider2 + '\n')
        
        # str_center = str_circle(ra_center, dec_center, 10)
        # f.write(str_center + '\n')
        
        str_ifs = str_box(ra_ifs, dec_ifs, x_ifs.value, y_ifs.value, pa_ifs+rot_angle)
        f.write(str_ifs + '\n')
        
        str_guider1 = str_box(ra_guider1, dec_guider1, x_guider1.value, y_guider1.value, pa_guider1+rot_angle)
        f.write(str_guider1 + '\n')

def str_box(ra, dec, x, y, angle):
    
    return 'box(' + ra + ',' + dec + ',' + str(x) + '",' + str(y) + '",' + str(angle) + ')'
    
def ifs2others(ra_ifs, dec_ifs, rot_angle, outfile):

    ra_gold, dec_gold = ifs2gold(ra_ifs, dec_ifs, rot_angle)
    
    ra_gnew, dec_gnew = gold2gnew(ra_gold, dec_gold, rot_angle)

    write_reg(outfile, ra_gnew, dec_gnew, ra_gold, dec_gold,
              ra_ifs, dec_ifs, rot_angle)


def ifs2gold(ra_ifs, dec_ifs, rot_angle):
    
    c_ifs = SkyCoord(ra_ifs, dec_ifs, unit=(u.hourangle, u.deg))

    position_angle = (rot_angle + pa_guider1) * u.deg
    separation = distance_guider1_ifs

    c_guider = c_ifs.directional_offset_by(position_angle, separation)  
    ra_guider = ra_string_format(c_guider.ra.to_string(u.hour))
    dec_guider = dec_string_format(c_guider.dec.to_string(u.degree))

    return ra_guider, dec_guider

def gold2gnew(ra_guider1, dec_guider1, rot_angle):
    
    c_guider1 = SkyCoord(ra_guider1, dec_guider1, unit=(u.hourangle, u.deg))

    position_angle = (271.07611342 + rot_angle) * u.deg + 180 * u.deg
    separation = 494.574 * u.arcsec

    c_guider2 = c_guider1.directional_offset_by(position_angle, separation)  
    
    ra_guider2 = ra_string_format(c_guider2.ra.to_string(u.hour))
    dec_guider2 = dec_string_format(c_guider2.dec.to_string(u.degree))

    return ra_guider2, dec_guider2



    
if __name__ == "__main__":

    
    rot_angle = 0   # 视场旋转角，单位为度。正北为0度，逆时针方向增加。   
    ra_ifs = '12:10:32.58'
    dec_ifs = '+39:24:21.1'
    ifs2others(ra_ifs, dec_ifs, rot_angle, outfile='NGC4151.reg')
    
    

        
    
    
