import Functions
import cProfile
import pstats

cProfile.run('Functions.currentDensity(0.1,0.1,1,0.5,500.,100.,0.,1,\'y1y2\')', 'profile.log')
p = pstats.Stats('profile.log')
p.sort_stats('tottime').print_stats(20)
p.sort_stats('cumtime').print_stats(20)