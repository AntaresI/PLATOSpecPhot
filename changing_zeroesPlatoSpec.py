    
import numpy as np

'''FUNCTION THAT CHANGES ZEROES TO AN AVERAGE OF SURROUNDING PIXELS IN AN IMAGE'''
def changing_zeroes(data):
        
        zeroinds = np.argwhere((data==0))   # zeroinds = np.argwhere((data==0) | (data==65535)) 
        zeroinds = tuple(map(tuple, zeroinds))       #convert array to tuple
        print("Indices where ADU is zero are: ", zeroinds)    
        if len(zeroinds) != 0:
             
                for i in range(len(zeroinds)):
                    value = 0 
                    counterval = 0
                    for j in range(zeroinds[i][0]-1, zeroinds[i][0]+2):
                        for k in range(zeroinds[i][1]-1, zeroinds[i][1]+2):
                            
                            if j >= data.shape[1]:
                                j -= 5
                            if k >= data.shape[1]:
                                k -= 5
                                
                            if data[j,k] != 0:
                                value += data[j,k]
                                counterval += 1
                                
                    value = int(value/counterval)
                    data[zeroinds[i]] = value 
        return data           
        ''''''


