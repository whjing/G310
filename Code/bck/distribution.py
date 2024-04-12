#%%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
import math
from mpl_toolkits.axisartist.axislines import AxesZero

def hii_circle(ll, bb, rad):
    radi = rad # in deg
    ang = np.arange(0.0,361.0,0.3)/180.0*3.1415926
    xx = radi*np.sin(ang)+ll
    yy = radi*np.cos(ang)+bb
    return xx,yy

def bar_pos(R0,ang,rad_c,bar_length):
    fac_0 = R0*np.cos(ang*math.pi/180.0)
    fac_1 = rad_c*rad_c - R0*R0 + fac_0**2.0
    dsource = np.sqrt(fac_1) + fac_0
    x_s = dsource*np.sin(ang*math.pi/180.0)
    y_s = R0-dsource*np.cos(ang*math.pi/180.0)

    dsou_c = dsource-bar_length
    
    x_sc = dsou_c*np.sin(ang*math.pi/180.0)
    y_sc = R0-dsou_c*np.cos(ang*math.pi/180.0)
    
#    print(x_s)
#    print(y_s)
#    print(dsource)
    return x_s,y_s,x_sc,y_sc


# ------- Data ---------
x = [-10,-5,0,5,10]
y = [3.0,-3.0,-2.0,-4.0,3.0]

file_factor2 = '../Data/coord-dis.txt'
with open(file_factor2,"r") as file:
    xx = len(file.readlines())
    lon = [0.000]*xx
    lat = [0.000]*xx
    dis = [0.000]*xx
    dise1 = [0.000]*xx
    dise2 = [0.000]*xx
print(xx)

with open(file_factor2,"r") as file2:
    kk = 0
    for line in file2:
        sou = line
        sou_a = sou.split('\n')
        sou_b = sou_a[0]
        sou_b0 = ' '.join(sou_b.split())
        sou_c = sou_b0.split(' ')

        lon[kk] = float(sou_c[0])
        lat[kk] = float(sou_c[1])
        dis[kk] = float(sou_c[2])
        dise1[kk] = float(sou_c[3])
        dise2[kk] = float(sou_c[4])

        kk = kk + 1
print('================')

#print(sdis)
#print(lon)



# -------- fig setup -----
# ------------------------
# ------------------------
spec_fig = '../Output/dis.png'


dirpath = '../Data/MW_2020.jpg'
img_plt = plt.imread(dirpath)
print("img_plt :",img_plt .shape)

#extent = (-24.4, 24.4, -24.4, 24.4)
extent = (-25.1, 25.1, -25.1, 25.1)  # much better
# ======================== box edge
x_lim0 = -15.5
x_lim1 = 15.5
y_lim0 = -15.5
y_lim1 = 15.5

# ======================== big circle
rad_plt = 15.0
c_linw = 1.0

# ======================== solor & GC
R0 = 8.5
solar_x = 0.0
solar_y = R0
GC_x = 0.0
GC_y = 0.0
s_size = 5
c_color = 'white'

txt_size = 7.0
txt_color = 'white'
txt_sun_x = -0.46
txt_sun_y = 7.5

txt_GC_x = 0.1
txt_GC_y = -0.7

# ======================== c-bar
bar_length = 0.5
bar_color = 'white'
bar_thick = 0.8
bar_fontsiz = 7.0 
bar_length_forc = bar_length+0.6

# ======================== axises acrossing the Sun
cro_color = 'gainsboro'
cro_thick = bar_thick*0.3
suncirc_thick = bar_thick*0.6

# ======================== HII setup
alpha_hii = 0.5
s_hii = 6.5
c_hii = 'white'
hii_color = 'white'
hii_thick = 0.8

s_hiio = 5.1
c_hiio = 'blue'

# ------------------------
# ------------------------
# -------- fig setup -----
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
fig, ax = plt.subplots(figsize=(5, 5), dpi=300.0)

plt.imshow(img_plt , cmap=plt.cm.binary, extent=extent)

ax.set_xlabel(r'$X (kpc)$', fontsize=15)
ax.set_ylabel(r'$Y (kpc)}$', fontsize=15)
#ax.set_title('Volume and percent change')

c_cw = 'white'
c_cw_overlay = 'fuchsia'
c_ccwb = 'black'
c_ccw = 'yellow'
marker_wb = 1.45
marker_w = 1.1
s_size_cal = 3.7
pii = 180.0/3.1415926

s_sh = 1.0
for i in range(len(lon)):
    Ds = dis[i]
    el = lon[i]
    be = lat[i]
    Ds_e1 = dise1[i]
    Ds_e2 = dise2[i]
    if (Ds > 0.0 and np.abs(Ds_e1) <5.0 and np.abs(Ds_e2) <5.0 ):
        sou_x=Ds*np.sin(el/pii)*np.cos(be/pii)
        sou_y=R0-Ds*np.cos(el/pii)*np.cos(be/pii)
        
        ax.plot(sou_x, sou_y, markersize = s_sh*s_size_cal, c=c_ccwb,marker="o",fillstyle='none',markeredgewidth=marker_wb)
        ax.plot(sou_x, sou_y, markersize = s_sh*s_size_cal, c=c_ccw,marker="o",fillstyle='none',markeredgewidth=marker_w)

        D1 = Ds + Ds_e1
        sx1 = D1*np.sin(el/pii)*np.cos(be/pii)
        sy1 = R0-D1*np.cos(el/pii)*np.cos(be/pii)
        
        D2 = Ds + Ds_e2
        sx2 = D2*np.sin(el/pii)*np.cos(be/pii)
        sy2 = R0-D2*np.cos(el/pii)*np.cos(be/pii)

        ax.plot([sx1,sx2],[sy1,sy2],c=c_ccw,linewidth=0.7)
        
xxc,yyc=hii_circle(0.0, 0.0, rad_plt)
ax.plot(xxc, yyc, c='white', linewidth=c_linw, linestyle='solid')

ax.set_xlim(x_lim0,x_lim1)
ax.set_ylim(y_lim0,y_lim1)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)
color_axis = 'white'
ax.spines['top'].set_color(color_axis)
ax.spines['left'].set_color(color_axis)
ax.spines['right'].set_color(color_axis)
ax.spines['bottom'].set_color(color_axis)


#===================================== cross the Sun axses:
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
txt_co = 'silver'

#xs,ys,xs_cr,ys_cr = bar_pos(R0,90.0,rad_plt,bar_length)
#xs,ys,xs_cl,ys_cl = bar_pos(R0,270.0,rad_plt,bar_length)
xx = [-rad_plt,rad_plt]
yy = [0.0,0.0]
ax.plot(xx,yy,color=txt_co,linewidth=cro_thick,linestyle='solid')#(0,(5,1)))

xs,ys,xs_cr,ys_cr = bar_pos(R0,0.0,rad_plt,bar_length)
xs,ys,xs_cl,ys_cl = bar_pos(R0,180.0,rad_plt,bar_length)
xx = [xs_cr,xs_cl]
yy = [ys_cr,ys_cl]
ax.plot(xx,yy,color=txt_co,linewidth=cro_thick,linestyle='solid')#(0,(5,1)))

#====== label on axises:
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
y_len = 0.15
s_ylen = 2.2

xx = [-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
#[-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6,7,8,9,10,11,12]
for i in range(len(xx)):
    xx0=[xx[i],xx[i]]
    yy0=[-y_len,y_len]
    if (xx[i] % 5 ==0):
        yy0=[-y_len*s_ylen,y_len*s_ylen]
    ax.plot(xx0,yy0,color=cro_color,linewidth=cro_thick,linestyle='solid')

#yy = [-23,-22,-21,-20,-19,-18,-17,-16,-15,-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,1,2,3,4,5,6]
yy = [-14,-13,-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14]
for i in range(len(yy)):
    xx0=[-y_len,y_len]
    if (yy[i] % 5 ==0):
        xx0=[-y_len*s_ylen,y_len*s_ylen]
    yy0=[yy[i],yy[i]]
    ax.plot(xx0,yy0,color=cro_color,linewidth=cro_thick,linestyle='solid')

y_shift = 0.9
txt_si = 5.7


txt_y = y_shift-1.7
txt_x = 5.0
ax.text(txt_x,txt_y, '5', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)
txt_x = 10.0
ax.text(txt_x,txt_y, '10', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)
txt_x = 10.9
ax.text(txt_x,txt_y, 'kpc', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)
txt_x = -5.0
ax.text(txt_x,txt_y, '5', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)
txt_x = -10.0
ax.text(txt_x,txt_y, '10', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)

txt_x = y_shift-0.3
txt_y = 5.0-0.2
ax.text(txt_x,txt_y, '5', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)
txt_y = -5.0-0.2
ax.text(txt_x,txt_y, '5', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)
txt_y = -10.0-0.2
ax.text(txt_x,txt_y, '10', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)

txt_x = 1.5
ax.text(txt_x,txt_y, 'kpc', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)

txt_x = y_shift-0.3
txt_y = 10.0-0.2
ax.text(txt_x,txt_y, '10', fontsize=txt_si,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='center',color=txt_co)



#======  Solar circle
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
x_sun,y_sun = hii_circle(0.0, 0.0, R0)
ax.plot(x_sun,y_sun,color='white',linewidth=suncirc_thick,linestyle=(0,(5,5)))#'dashed')


#========================================================================= bar-plotting
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
ang_list=[10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,120.0,150.0]
# ===== quadrant-1
for i in range(len(ang_list)):
    ang = ang_list[i]
    xs,ys,xs_c,ys_c = bar_pos(R0,ang,rad_plt,bar_length)
    xx = [xs_c,xs]
    yy = [ys_c,ys]
    ax.plot(xx,yy,color=bar_color,linewidth=bar_thick)
    ang_txt = format(ang,'.0f')+r"$^{\circ}$"
    xs,ys,xs_c,ys_c = bar_pos(R0,ang,rad_plt,bar_length_forc)
    ax.text(xs_c,ys_c, ang_txt, fontsize=bar_fontsiz,
              rotation=90.0+ang+180, rotation_mode='anchor',horizontalalignment='center',color=bar_color)

ang_list=[210.0,240.0,270.0,280.0,290.0,300.0,310.0,320.0,330.0,340.0,350.0]
# ===== quadrant-4
for i in range(len(ang_list)):
    ang = ang_list[i]
    xs,ys,xs_c,ys_c = bar_pos(R0,ang,rad_plt,bar_length)
    xx = [xs_c,xs]
    yy = [ys_c,ys]
    ax.plot(xx,yy,color=bar_color,linewidth=bar_thick)
    ang_txt = format(ang,'.0f')+r"$^{\circ}$"
    xs,ys,xs_c,ys_c = bar_pos(R0,ang,rad_plt,bar_length_forc+0.2)
    ax.text(xs_c,ys_c, ang_txt, fontsize=bar_fontsiz,
              rotation=90.0+ang, rotation_mode='anchor',horizontalalignment='center',color=bar_color)
    
ang_list=[0.0,180.0]
for i in range(len(ang_list)):
    ang = ang_list[i]
    xs,ys,xs_c,ys_c = bar_pos(R0,ang,rad_plt,bar_length)
    xx = [xs_c,xs]
    yy = [ys_c,ys]
    ax.plot(xx,yy,color=bar_color,linewidth=bar_thick)


#===================================== Sun & GC
# - -           ----------------------
# - -           ----------------------
# - -           ----------------------
ax.plot(solar_x, solar_y, markersize = s_size, c=c_color, marker=(5, 1))
ax.plot(solar_x, solar_y, markersize = s_size*0.8, c='black', marker=(5, 1))

ax.plot(GC_x, GC_y, markersize = s_size,c=c_color, marker=(5, 1))
ax.plot(GC_x, GC_y, markersize = s_size*0.8,c='black', marker=(5, 1))

ax.text(txt_sun_x,txt_sun_y, 'Sun', fontsize=txt_size,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='left',color=txt_color)

ax.text(txt_GC_x,txt_GC_y, 'GC', fontsize=txt_size,
              rotation=0.0, rotation_mode='anchor',horizontalalignment='left',color=txt_color)   

#ax.scatter(GC_x, GC_y, s=s_size, c=c_color, marker=(5, 1))



#ax.grid(True)
fig.tight_layout()


ax.set_box_aspect(1)
plt.savefig(spec_fig, bbox_inches='tight',format="png")
plt.close('all')
#plt.show()

# %%
