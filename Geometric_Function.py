# -*- coding: utf-8 -*-
def convertPolygonPtLsttoLn(ptLst):
    """
    ptLst=[[45027.4,-76720.2],
           [-5522.6,-76720.2],
           [-5522.6,-9549.1],
           [3538.0,-1100.0],
           [45027.4,-1100.0]]    
    """
    lnLst=[]
    for i in range(len(ptLst)-1):
        lnLst.append([ptLst[i],ptLst[i+1]])
    lnLst.append([ptLst[-1],ptLst[0]])
    return lnLst

def ptInPolygon(pt,ptLst):
    """
    pt=[-6000,-77000]
    ptLst=[[45027.4,-76720.2],
           [-5522.6,-76720.2],
           [-5522.6,-9549.1],
           [3538.0,-1100.0],
           [45027.4,-1100.0]]
    or 
    ptLst=((45027.4,-76720.2),
           (-5522.6,-76720.2),
           (-5522.6,-9549.1),
           (3538.0,-1100.0),
           (45027.4,-1100.0))
    """
    if(type(ptLst)==tuple):
        ptLst=list(ptLst)
        for i in range(len(ptLst)):
            ptLst[i]=list(ptLst[i])
    inPolygon=False
    xPt=pt[0]
    yPt=pt[1]
    lnLst=convertPolygonPtLsttoLn(ptLst)
    xExt=[]
    yExt=[]
    for i in range(len(lnLst)):
        xExt.append((min(lnLst[i][0][0],lnLst[i][1][0]),max(lnLst[i][0][0],lnLst[i][1][0])))
        yExt.append((min(lnLst[i][0][1],lnLst[i][1][1]),max(lnLst[i][0][1],lnLst[i][1][1])))
    intCount=0
    for i in range(len(lnLst)):
        #if(yPt>=yExt[i][0] and yPt<yExt[i][1] and xPt<=xExt[i][1]):
        if(((yPt>=yExt[i][0] and yPt<yExt[i][1]) or (yPt==yExt[i][0] and yPt==yExt[i][1])) and xPt<=xExt[i][1]): 
            xPtInterp=interpForX(lnLst[i],yPt)
            if(xPt<xPtInterp or pt in ptLst):
                intCount+=1
    if(intCount>0 and intCount%2==1):
        inPolygon=True
    return inPolygon

def interpForX(ptList,yPt):
    """ptList is the 2 points that make up the line
    eg.ptList=[[-310,-67.5],[-310,-837.5]]
    """
    delX=ptList[1][0]-ptList[0][0]
    delY=ptList[1][1]-ptList[0][1]
    if(delY==0):
        delY=0.0000000001
    changeY=yPt-ptList[0][1]
    xPt=ptList[0][0]+changeY/delY*delX
    return xPt

def polygonArea(ptList):
    '''
    ptList is a list of points example
    ptList=[[25,0],
            [775,0],
            [800,25],
            [800,225],
            [775,250],
            [25,250],
            [0,225],
            [0,25]]
    '''
    areaSum=0
    centroidXSum=0
    centroidYSum=0
    for i in range(len(ptList)):
        if(i==(len(ptList)-1)):
            area=ptList[i][0]*ptList[0][1]-ptList[0][0]*ptList[i][1]
            centroidX=(ptList[i][0]+ptList[0][0])*area
            centroidY=(ptList[i][1]+ptList[0][1])*area
        else:
            area=ptList[i][0]*ptList[i+1][1]-ptList[i+1][0]*ptList[i][1]
            centroidX=(ptList[i][0]+ptList[i+1][0])*area
            centroidY=(ptList[i][1]+ptList[i+1][1])*area
        areaSum=areaSum+area
        centroidXSum+=centroidX
        centroidYSum+=centroidY
    pgArea=0.5*abs(areaSum)
    centroidX=centroidXSum/(6*pgArea)
    centroidY=centroidYSum/(6*pgArea)
    return pgArea,centroidX,centroidY