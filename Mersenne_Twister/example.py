import Mersenne_Twister

MT = Mersenne_Twister.Mersenne_Twister()# if no args , seeds length is 4 and each seed number is initialized by time
print MT.get_32(),MT.get_31(),MT.get_real1(),MT.get_real2(),MT.get_real3()
MT1 = Mersenne_Twister.Mersenne_Twister(5)# one args(int) is seeds length , each seed number is initialized by time
print MT1.get_32(),MT1.get_31(),MT1.get_real1(),MT1.get_real2(),MT1.get_real3()
MT2 = Mersenne_Twister.Mersenne_Twister([0x123, 0x234, 0x345, 0x456]) #one args(list) is seeds
print MT2.get_32(),MT2.get_31(),MT2.get_real1(),MT2.get_real2(),MT2.get_real3()

