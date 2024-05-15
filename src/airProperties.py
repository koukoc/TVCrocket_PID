import numpy as np
class AirProperties:
    def __init__(self,height):
        assert type(height) is float or int,'wrong input'
        self.height = height
        self.getProperties()

    def getProperties(self):
        if self.height < 11000:
            self.Temperature = 15.04-0.00649*self.height
            self.pressure = 101.29 * np.power(((self.Temperature + 273.1)/288.08),5.256)
        
        if self.height >= 11000 and self.height < 25000:
            self.Temperature = -56.46
            self.pressure = 22.65 * np.exp(1.73-0.000157*self.height)
        
        if self.height >= 25000:
            self.Temperature = -131.21+0.00299*self.height
            self.pressure = 2.488 * np.power(((self.Temperature + 273.1)/216.6),-11.388)
        
        self.density = self.pressure/(0.2869*(self.Temperature+273.1))