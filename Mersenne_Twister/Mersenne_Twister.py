"""
A Python-program for MT19937

I reffered to "http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/mt19937ar.html"

This code License is BSD in accordance with "http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/license.html"


"""

class Mersenne_Twister():
    __N = 624 #default
    __M = 397 #default
    __genrand_number = 19650218 #default
    __MATRIX_A  = 0x9908b0df
    __mag01 = [0,__MATRIX_A]
    __UPPER_MASK = 0x80000000
    __LOWER_MASK = 0x7fffffff
    def __init__(self,*seeds):
        if len(seeds) == 1 and type(seeds[0]) == list:
            init_list = seeds[0]
        elif (len(seeds) == 1 and type(seeds[0]) == int) or (len(seeds) == 0):
            import time
            if len(seeds) == 0:
                init_list_length = 4
            else:
                init_list_length = seeds[0]
            init_list = [0 for _ in xrange(init_list_length)]
            for i in range(init_list_length):
                init_list[i] = int(time.time() * (10**2) + i) & 0xffffffff
        else:
            raise TypeError,"Mersenne twister init by list"

        self._init_genrand(self.__genrand_number)

        i = 1
        j = 0
        k = max(self.__N,len(init_list))
        for _ in xrange(k):
            self.mt[i] = (self.mt[i] ^ (self.mt[i-1] ^ (self.mt[i-1] >> 30)) * 1664525) + init_list[j] + j
            self.mt[i] &= 0xffffffff
            i += 1
            j += 1
            if i >= self.__N:
                self.mt[0] = self.mt[self.__N-1]
                i = 1
            if j >= len(init_list):
                j = 0
        for _ in xrange(k-1):
            self.mt[i] = (self.mt[i] ^ ((self.mt[i-1] ^ (self.mt[i-1] >> 30)) * 1566083941)) - i
            self.mt[i] &= 0xffffffff
            i += 1
            if i >= self.__N:
                self.mt[0] = self.mt[self.__N-1]
                i = 1
        self.mt[0] = 0x80000000
        
    def _init_genrand(self,init_number):
        self.mt = [0 for _ in xrange(self.__N)]
        self.mt[0] = init_number & 0xffffffff
        for i in range(1,self.__N):
            self.mt[i] = 1812433253 * (self.mt[i-1] ^ (self.mt[i-1] >> 30)) + i
            self.mt[i] &= 0xffffffff
        self.mti = self.__N
            
    def get_32(self):#return [0,0xffffffff] range unsigned int 
        if self.mti >= self.__N:
            
            if self.mti == self.__N+1:#init_genrand has not been called
                self._init_genrand(5489)
            
            for i in range(self.__N-self.__M):
                y = (self.mt[i]&self.__UPPER_MASK) | (self.mt[i+1]&self.__LOWER_MASK)
                self.mt[i] = self.mt[i+self.__M] ^ (y >> 1) ^ self.__mag01[y & 0x1]
            
            for i in range(self.__N-self.__M,self.__N-1):
                y = (self.mt[i]&self.__UPPER_MASK) | (self.mt[i+1]&self.__LOWER_MASK)
                self.mt[i] = self.mt[i+(self.__M-self.__N)] ^ (y >> 1) ^ self.__mag01[y&0x1]
            
            self.mti = 0
        
        y = self.mt[self.mti]
        self.mti += 1
        
        y ^= (y >> 11)
        y ^= (y << 7) & 0x9d2c5680
        y ^= (y << 15) & 0xefc60000
        y ^= (y >> 18)
        return y

    def get_31(self): #return [0,0x7fffffff] range signed int
        return (self.get_32() >> 1)
        
    def get_real1(self): #return [0,1] real number
        return (self.get_32() * (1.0/4294967295.0))
    def get_real2(self): #return [0,1) real number
        return (self.get_32() * (1.0/4294967296.0))
    def get_real3(self): #return (0,1) real number
        return (self.get_32() + 0.5) * (1.0/4294967296.0)
