from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.metrics import mean_squared_error
from numpy import *
import numpy as np
from matplotlib import *

def normalizer(data):
    normalized_data = np.divide(data-np.mean(data),np.max(data)-np.min(data))
    return(normalized_data)
  
def Xsm(X_train,y_train):
    degrees = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
    N=len(degrees)
    P=np.empty(N)
    xsm=np.empty(N)
    pred_y=np.empty([N,len(y_train)])
    for i in range(N):
        # Calculate Xsm
        polynomial_features = PolynomialFeatures(degree=degrees[i])
        x_poly = polynomial_features.fit_transform(X_train)
        model_L = LinearRegression()
        model_L.fit(x_poly, y_train)
        y_poly_pred = model_L.predict(x_poly)
        xsm[i] = mean_squared_error(y_train,y_poly_pred)
    
    return(degrees,xsm)
        
def Zsm(xsm,degrees):
    M=len(degrees)
    VarDs =arange(0.1,1,0.01)
    alpha,beta=5,5
    final_zm=np.ones(shape=[M,len(VarDs)])*100
    
    s2=0
    
    for sigma in VarDs:
        mn = 0
        s1=0
        for m in degrees:
            mw=(1-np.divide(m,M))*(sigma**2)
            Km=np.divide((2*alpha*sigma),np.sqrt(M))*(np.sqrt((np.divide((alpha**2)*(sigma**2),M)+xsm[mn]-(0.5*mw))))
            zsmu=xsm[mn]-mw+(np.divide((2*(alpha**2)*(sigma**2)),M))+Km+(np.divide(m,M)*sigma**2)+beta*(np.sqrt(2*m)/M)*sigma**2
            mn = mn+1
            if np.iscomplex(zsmu):
                pass
            else:
                final_zm[s1,s2] = zsmu
                s1=s1+1
        s2=s2+1
        
    ind1,ind2 = np.unravel_index(np.argmin(final_zm, axis=None), final_zm.shape)
    order_opt = degrees[ind1] 
    var_opt = VarDs[ind2]
        
   # print("Optimum varaince: "+str(var_opt))    
           
    return(order_opt,final_zm)

def main(X_train,y_train):
    X_train = normalizer(X_train)
    y_train = normalizer(y_train)
    
    degrees,xsm_val=Xsm(X_train,y_train)
    order_opt,final_zm=Zsm(xsm_val,degrees)
    return(order_opt,final_zm)
