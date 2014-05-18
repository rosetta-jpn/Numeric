import scipy.linarg 
import numpy as np
class RangeError(Exception):
    pass
    
class Poisson():
    __h = 0.5
    __EPS = 10 ** (-8)
    def __init__(self,rangex,rangey,f,Dirichlet,Neuman):
        self.f = f
        self.Dirichlet = Dirichlet
        self.Neuman = Neuman
        if self._rangeOK(rangex):
            self.rangex = rangex
        else:
            raise RangeError,"rangex"
        if self._rangeOK(rangey):
            self.rangey = rangey
        else:
            raise RangeError,"rangey"

        self._solver()
        self.xmax = int((self.rangex[1]-self.rangex[0]) / self.__h)
        self.ymax = int((self.rangey[1]-self.rangey[0]) / self.__h)

    def _rangeOK(self,r):
        return type(rangex) == list and len(rangex) != 2 and type(rangex[0]) == type(rangex[1]) and type(rangex[0]) == int and rangex[0] < rangex[1]
        
    def calc(self,x,y):
        if x % self.__h < self.__EPS or (self.__h - (x % self.__h)) < self.__EPS:
            _x = (x - self.rangex[0]) / self.__h
        if y % self.__h < self.__EPS or (self.__h - (y % self.__h)) < self.__EPS:
            _y = (y-self.rangex[0]) / self.__h
        
        if not(0 <= _x <= self.xmax and 0 <= _y <= self.ymax):
            print "(%f,%f) is not calced" % (x,y)
        else:
            print "u(%f,%f) = %f" % (x,y,self.ans[_x][_y])

            
    def _solver(self):
        #scipy
        Matrix_w = int ((self.rangey[1]-self.rangey[0]) /self.__h)
        Matrix_h = int ((self.rangex[1]-self.rangex[0]) /self.__h)
        A = np.zeros((Matrix_h,Matrix_w))
        for i in range(Matrix_h):
            for j in range(Matrix_w):
                pass
        b = np.zeros(Matrix_w)
        for i in range(Matrix_w):
            pass
        
        x = scipy.linarg.solve(A,b)
        self.ans = [[0 for _ in xrange(Matrix_w)] for _ in xrange(Matrix_h)]
        for i in range(Matrix_h):
            for j in range(Matrix_w):
                ans[i][j] = x[i*Matrix_h+j]
        
