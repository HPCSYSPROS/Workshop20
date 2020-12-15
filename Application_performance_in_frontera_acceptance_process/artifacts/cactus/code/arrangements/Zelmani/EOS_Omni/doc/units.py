# Copyright Christian D. Ott, Dec 29, 2012
# cott@tapir.caltech.edu
#
# python script that takes fundamental constants and spits out
# conversion factors between c = G = M_Sun = 1 and cgs

# constants from NIST table, 2010 data
# solar mass from http://asa.usno.navy.mil/SecK/2013/Astronomical_Constants_2013.pdf

ggrav = 6.6738480e-8 #cgs
clite = 2.99792458e10 #cgs
msun = 1.9884e33 #cgs

print "Using the following constants:"
print "G: 6.6738480e-8 cm**3 / g / s**2"
print "c: 2.99792458e10 cm / s"
print "M_sun: 1.9884e33 g"

print "      mass_gf=%22.14E" % (1.0/msun)
print "  inv_mass_gf=%22.14E" % msun

# length
length_gf =  1.0/(ggrav/clite**2 * msun)
inv_length_gf = 1.0/length_gf

print "    length_gf=%22.14E" % length_gf
print "inv_length_gf=%22.14E" % inv_length_gf

# density
rho_gf = 1.0/msun * 1.0/length_gf**3
inv_rho_gf = 1.0/rho_gf

print "       rho_gf=%22.14E" % rho_gf
print "   inv_rho_gf=%22.14E" % inv_rho_gf

# time
time_gf = clite * length_gf
inv_time_gf = 1.0/time_gf
print "      time_gf=%22.14E" % time_gf
print "  inv_time_gf=%22.14E" % inv_time_gf

# specific internal energy erg / gram
eps_gf = 1.0/clite**2
inv_eps_gf = clite**2
print "       eps_gf=%22.14E" % eps_gf
print "   inv_eps_gf=%22.14E" % inv_eps_gf

# pressure
press_gf = 1.0/msun * inv_time_gf**2 * inv_length_gf
inv_press_gf = 1.0/press_gf
print "     press_gf=%22.14E" % press_gf
print " inv_press_gf=%22.14E" % inv_press_gf

