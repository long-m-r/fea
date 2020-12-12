import logging
log=logging.getLogger('fea.util')

import numpy as np

#>>>>>>>>>>>>>>>>>> UTILITIES <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
def lstsq(A,b,eps,title="LstSq Solve"):
    
    try:
        # If square A, use np.linalg.solve instead of np.linalg.lstsq. This should be faster (I hope!)
        shape=np.shape(A)
        if shape[0]==shape[1]:
            res=np.linalg.solve(A,b)
            if log.getEffectiveLevel()<=10: # VERY expensive logging operation, but incredibly useful
                _log_lstsq(A,b,[res,[]],title='Direct Solve')
            return res
        else:
            
            res=np.linalg.lstsq(A,b)
            if log.getEffectiveLevel()<=10: # VERY expensive logging operation, but incredibly useful
                _log_lstsq(A,b,res,title=title)


            if np.sum(np.power(np.dot(A,res[0])-b,2))>len(A)*eps**2:
                raise ValueError('LstSq Result was invalid')
            return res[0]
    except np.linalg.linalg.LinAlgError:
        raise ValueError('LstSq Result was invalid')
        

def _log_lstsq(A,b,res,pre='\t',title="LstSq Solve"):
    # log.info('Solving lstsq:')
    rstr=title+" (residuals: "+str(res[1])+"):\n"
    h,w=np.shape(A)

    for i in range(max(h,w)):
        mida=" "
        midb=" "
        astart="│"
        aend="│"
        xstart="|"
        xend="|"

        if i==0:
            midb="="
            mida="x"
            if h==1:
                astart="["
                aend="]"
            else:
                astart="┌"
                aend="┐"
            if w==1:
                xstart="["
                xend="]"
            else:
                xstart="┌"
                xend="┐"
        else:
            if i==h-1:
                astart="└"
                aend="┘"
            elif i>=h:
                astart=" "
                aend=" "

            if i==w-1:
                xstart="└"
                xend="┘"
            elif i>=w:
                xstart=" "
                xend=" "

        if i>=w:
            rstr+=pre+astart+",".join(["{:8.2f}".format(a) for a in A[i]])+aend+mida+xstart+"        "+xend+midb+astart+"{:8.2f}".format(b[i])+aend+' ('+"{:8.2f}".format((np.dot(res[0],A[i])-b[i])**2)+')\n'
        elif i>=h:
            rstr+=pre+astart+" ".join(["        ".format(a) for a in A[0]])+aend+mida+xstart+"{:8.2f}".format(res[0][i])+xend+midb+astart+"        "+aend+'\n'
        else:
            rstr+=pre+astart+",".join(["{:8.2f}".format(a) for a in A[i]])+aend+mida+xstart+"{:8.2f}".format(res[0][i])+xend+midb+astart+"{:8.2f}".format(b[i])+aend+' ('+"{:8.2f}".format((np.dot(res[0],A[i])-b[i])**2)+')\n'
    log.info(rstr[:-1])
