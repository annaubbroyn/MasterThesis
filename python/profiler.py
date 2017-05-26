import Functions
import cProfile
import pstats

cProfile.run('Functions.dFreeEnergy(0.5,0.0,0.2,0.1,500.,100.,0,1.,\'y1y2\')', 'profile.log')
p = pstats.Stats('profile.log')
p.sort_stats('cumtime').print_stats(20)