from astropy import units as u
from astropy.coordinates import SkyCoord
import pdb

x_guider2 = 648 * u.arcsec
y_guider2 = 442.2 * u.arcsec
angle_guider2 = 320

x_guider1 = 318 * u.arcsec
y_guider1 = 240 * u.arcsec

x_ifs = 65 * u.arcsec       # -2 degree ~ 22.6 arcsec @ 10.8 arcmin
y_ifs = 71 * u.arcsec
 

def gold2ifs(ra_guider, dec_guider):

    c_guider = SkyCoord(ra_guider, dec_guider, unit=(u.hourangle, u.deg))

    position_angle = -2 * u.deg + 180 * u.deg
    separation = 651.563 * u.arcsec

    c_ifs = c_guider.directional_offset_by(position_angle, separation)  

    ra_ifs = ra_string_format(c_ifs.ra.to_string(u.hour))
    dec_ifs = dec_string_format(c_ifs.dec.to_string(u.degree))

    return ra_ifs, dec_ifs


def ifs2gold(ra_ifs, dec_ifs):
    
    c_ifs = SkyCoord(ra_ifs, dec_ifs, unit=(u.hourangle, u.deg))

    position_angle = -2 * u.deg
    separation = 651.563 * u.arcsec

    c_guider = c_ifs.directional_offset_by(position_angle, separation)  
    ra_guider = ra_string_format(c_guider.ra.to_string(u.hour))
    dec_guider = dec_string_format(c_guider.dec.to_string(u.degree))

    return ra_guider, dec_guider



def gnew2center(ra_guider2, dec_guider2):
    
    c_guider2 = SkyCoord(ra_guider2, dec_guider2, unit=(u.hourangle, u.deg))

    position_angle = 50 * u.deg + 180 * u.deg
    separation = 766.5 * u.arcsec

    c_center = c_guider2.directional_offset_by(position_angle, separation)  
    
    ra_center = ra_string_format(c_center.ra.to_string(u.hour))
    dec_center = dec_string_format(c_center.dec.to_string(u.degree))

    return ra_center, dec_center

def center2ifs(ra_center, dec_center):
    
    c_center = SkyCoord(ra_center, dec_center, unit=(u.hourangle, u.deg))

    position_angle = 180 * u.deg 
    separation = 64.458 * u.arcsec

    c_ifs = c_center.directional_offset_by(position_angle, separation)  
    
    ra_ifs = ra_string_format(c_ifs.ra.to_string(u.hour))
    dec_ifs = dec_string_format(c_ifs.dec.to_string(u.degree))

    return ra_ifs, dec_ifs
    
    
def center2guider1(ra_center, dec_center):
    
    c_center = SkyCoord(ra_center, dec_center, unit=(u.hourangle, u.deg))

    position_angle = 0 * u.deg 
    separation = 651.563 * u.arcsec

    c_guider1 = c_center.directional_offset_by(position_angle, separation)  
    
    ra_guider1 = ra_string_format(c_guider1.ra.to_string(u.hour))
    dec_guider1 = dec_string_format(c_guider1.dec.to_string(u.degree))

    return ra_guider1, dec_guider1

def gnew2gold(ra_guider2, dec_guider2):
    
    c_guider2 = SkyCoord(ra_guider2, dec_guider2, unit=(u.hourangle, u.deg))

    position_angle = 271.07611342 * u.deg 
    separation = 478.669068 * u.arcsec

    c_guider1 = c_guider2.directional_offset_by(position_angle, separation)  
    
    ra_guider1 = ra_string_format(c_guider1.ra.to_string(u.hour))
    dec_guider1 = dec_string_format(c_guider1.dec.to_string(u.degree))

    return ra_guider1, dec_guider1


def gold2gnew(ra_guider1, dec_guider1):
    
    c_guider1 = SkyCoord(ra_guider1, dec_guider1, unit=(u.hourangle, u.deg))

    position_angle = 271.07611342 * u.deg + 180 * u.deg
    separation = 478.669068 * u.arcsec

    c_guider2 = c_guider1.directional_offset_by(position_angle, separation)  
    
    ra_guider2 = ra_string_format(c_guider2.ra.to_string(u.hour))
    dec_guider2 = dec_string_format(c_guider2.dec.to_string(u.degree))

    return ra_guider2, dec_guider2



def offsets(ra1, dec1, ra2, dec2):
    
    pos1 = SkyCoord(ra1, dec1, frame='icrs', unit=(u.hourangle, u.deg))
    pos2 = SkyCoord(ra2, dec2, frame='icrs', unit=(u.hourangle, u.deg))
    dra, ddec = pos1.spherical_offsets_to(pos2)
    print('delta_ra: ', dra.to(u.arcsec))
    print('delta_dec: ', ddec.to(u.arcsec))
    

def offsets2(ra1, dec1, ra2, dec2):
    
    ra1 = '13 42 44.0210'
    dec1 = '+28 26 18.340'
    
    ra2 = '13 42 07.7090'
    dec2 = '+28 26 26.388'
    
    pos1 = SkyCoord(ra1, dec1, frame='icrs', unit=(u.hourangle, u.deg))
    pos2 = SkyCoord(ra2, dec2, frame='icrs', unit=(u.hourangle, u.deg))
    sep = pos1.separation(pos2)
    pos = pos1.position_angle(pos2).to(u.deg)

def offsets3():
    
    ra1 = '15 52 18.4479'
    dec1 = '+33 01 33.828'
    dra = 220.964 * u.arcsec 
    ddec = 363.868 * u.arcsec
    
    pos1 = SkyCoord(ra1, dec1, frame='icrs', unit=(u.hourangle, u.deg))
    ra_off = pos1.ra + dra
    dec_off = pos1.dec + ddec 
    pos2 = SkyCoord(ra_off, dec_off, frame='icrs', unit=(u.deg, u.deg))
    print(ra_string_format(pos2.ra.to_string(u.hour)))
    print(dec_string_format(pos2.dec.to_string(u.degree)))

def str_box(ra, dec, x, y, angle):
    
    return 'box(' + ra + ',' + dec + ',' + str(x) + '",' + str(y) + '",' + str(angle) + ')'
    
    
def str_circle(ra, dec, radius):
    
    return 'circle(' + ra + ',' + dec + ',' + str(radius) + '")'


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


def write_reg(outfile, 
              ra_guider2, dec_guider2, 
              ra_ifs, dec_ifs,
              ra_guider1, dec_guider1):
    
    with open(outfile, "w") as f:
        f.write('# Region file format: DS9 version 4.1 \n')
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        f.write('fk5 \n')
        
        str_guider2 = str_box(ra_guider2, dec_guider2, 
                              x_guider2.value, y_guider2.value, 
                              angle=angle_guider2)
        f.write(str_guider2 + '\n')
        
        # str_center = str_circle(ra_center, dec_center, 10)
        # f.write(str_center + '\n')
        
        str_ifs = str_box(ra_ifs, dec_ifs, x_ifs.value, y_ifs.value, -2)
        f.write(str_ifs + '\n')
        
        str_guider1 = str_box(ra_guider1, dec_guider1, x_guider1.value, y_guider1.value, -2)
        f.write(str_guider1 + '\n')


def gnew2others(ra_gnew, dec_gnew, outfile):

    ra_guider1, dec_guider1 = gnew2gold(ra_gnew, dec_gnew)
    
    ra_ifs, dec_ifs = gold2ifs(ra_guider1, dec_guider1)

    write_reg(outfile, ra_gnew, dec_gnew, 
              ra_ifs, dec_ifs,
              ra_guider1, dec_guider1)
    

def ifs2others(ra_ifs, dec_ifs, outfile):

    ra_gold, dec_gold = ifs2gold(ra_ifs, dec_ifs)
    
    ra_gnew, dec_gnew = gold2gnew(ra_gold, dec_gold)

    write_reg(outfile, ra_gnew, dec_gnew, 
              ra_ifs, dec_ifs,
              ra_gold, dec_gold)
    

def offset4(ra_gnew, dec_gnew, ra_target, dec_target):
    
    ra_guider1, dec_guider1 = gnew2gold(ra_gnew, dec_gnew)
    
    ra_ifs, dec_ifs = gold2ifs(ra_guider1, dec_guider1)
    
    offsets(ra_ifs, dec_ifs, ra_target, dec_target)

    
if __name__ == "__main__":
    
    
    # ra_gnew = '15:52:18.4479'
    # dec_gnew = '+33:01:33.828'
    # gnew2others(ra_gnew, dec_gnew, 'bd32_test1.reg')


    ra_ifs = '10:39:36.7644'
    dec_ifs = '+43:06:08.896'
    ifs2others(ra_ifs, dec_ifs, 'feige34.reg')
    

    ra_gnew = '10:40:05.5619'
    dec_gnew = '+43:15:04.255'
    ra_target = '10:39:36.7644'
    dec_target = '+43:06:08.896'
    offset4(ra_gnew, dec_gnew, ra_target, dec_target)
    
    



