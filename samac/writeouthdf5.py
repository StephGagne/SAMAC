
import numpy as np
import os

#########################################################################        
###########################   writeouthdf5   ############################
#########################################################################
def writeouthdf5(self,fName=None, fdir=None):
    """This method writes out the content of the cloud in a string representation. This method can be used if one wants to export the cloud's data to a different programming language or make a language-independent backup of their cloud-object's data. Uses Python's repr function.
    fName: Name of the file where the cloud will be saved. Default saves with name SAMAC_Cloud_date_time.
    fdir: directory in which the file will be saved.
    WARNING: this method is designed to handle the current official structure of the object, if the structure has been customized by users, these custom structures will not be saved. This method will need to be modified to accomodate the new structure members.
    WARNING: any masked data will lost their masks. The feature to save masks is in planning for future versions."""
    #MMM writeout the masks of the arrays also!
    import Tkinter, tkFileDialog
    import h5py
    
    if fName==None: 
        if "layered" in self.desc:
            if self.desc["layered"]=='no': layno=''
            else: layno='_'+str(self.desc["layered"])
        else: layno=''
        fName="SAMAC_Cloud_" + self.desc["date"] + "_" + self.desc["time"] + layno      # name of the file where the cloud will be saved
    if fdir==None: fdir=os.getcwd()+'/'
    try:        # verifying if the file already exists
        f=open(fdir+fName+'.hdf5','r'); f.close();
        X=True
    except: X=False;
    Abort=False     # not aborting by default
    if X:
        print("File %s already exists in this directory (%s)." % (fName+'.hdf5',fdir))
        while X==True:
            Q=raw_input("Do you want to change the file name (N), change directory (D), both (ND) or overwrite (O)?  ")
            if 'n' in Q.lower():
                fName=str(raw_input("Enter the new file name:  "))
                try: 
                    f=open(fName+'.hdf5','r'); f.close();
                    X=True
                    print("File %s already exists in this directory (%s)." % (fName+'.hdf5',fdir))
                except: X=False
            if 'd' in Q.lower():
                root = Tkinter.Tk()
                print("Please choose the directory in which the file should be saved:")
                pathdef=""
                while len(pathdef)==0:
                    pathdef = tkFileDialog.askdirectory(parent=root,title='Please select a directory')
                    root.withdraw()
                    if len(pathdef) > 0:
                        print("You chose %s" % pathdef) 
                        fdir=pathdef+'/'
                        try:        # verifying if the file already exists
                            f=open(fdir+fName+'.hdf5','r'); f.close();
                            print("File %s already exists in this directory (%s)." % (fName+'.hdf5',fdir));
                            X=True
                        except: X=False;
                        break
                    else:
                        pd=raw_input("You must choose a directory. Would you like to choose again?  ")
                        if pd=='y': pass
                        else: break
            if Q.lower()=='o':
                X=False
            if Q.lower()=='o' or 'n' in Q.lower() or 'd' in Q.lower(): pass
            else: X=False; Abort=True
    if X: print("You are caught in a bug, try fixing the code or try something else. Sorry for the inconvenience.")
    else: 
        if Abort: print("Aborted.")
        else:
            # saving data in file
            f=h5py.File(fdir+fName+'.hdf5','w')
            data=f.create_dataset("data",data=self.data)
            desc=f.create_group("desc"); 
            for key in self.desc: 
                desc.attrs[key]=np.string_(self.desc[key])
            descedit=f.create_dataset("descedit",data=self.descedit)
            dttl=f.create_dataset("dttl",data=np.string_(self.dttl))
            dunit=f.create_dataset("dunit",data=np.string_(self.dunit))
            for i,actdata in enumerate(self.extradata):
                extradata=f.create_dataset("extradata"+str(i),data=actdata)
            xdnum=f.create_dataset("extradatanum",data=self.extradatanum)
            for i, actttl in enumerate(self.extrattl):
                xttl=f.create_dataset("extrattl"+str(i),data=np.string_(actttl))
            for i, actunit in enumerate(self.extraunit):
                xunit=f.create_dataset("extraunit"+str(i),data=np.string_(actunit))
            orient=f.create_dataset("orient",data=np.string_(self.orient))
            props=f.create_group("props");
            for key in self.props:
                if isinstance(self.props[key],np.ndarray):
                    dset=props.create_dataset(key,data=self.props[key])
                elif isinstance(self.props[key],list):
                    for i,L in enumerate(self.props[key]):
                        if isinstance(L,np.ndarray): dset=props.create_dataset(key+str(i),data=L)
                        elif isinstance(L,basestring): props.attrs[key+str(i)]=np.string_(L)
                        else: print("C.props[%s] is a list C.props[%s][%s] is neither an array nor a string" % (key,i))
                elif isinstance(self.props[key],basestring):
                    props.attrs[key]=np.string_(self.props[key])
                else: print("C.props[%s] not saved: type %s not handled" %(key,type(self.props[key])))
            for i,sd in enumerate(self.sd):
                szdst=f.create_group("sd"+str(i))
                for key in sd:
                    if isinstance(sd[key],np.ndarray):
                        dset=szdst.create_dataset(key,data=sd[key])
                    elif isinstance(sd[key],basestring):
                        szdst.attrs[key]=np.string_(sd[key])
                    else: print("C.sd[%s][%s] is of type %s, this type is not expected and not handled by this method." %(i,key,type(sd[key])))
            times=f.create_group("times")
            for key in self.times:
                dset=times.create_dataset(key,data=self.times[key])
            f.close()
            print("File saved at %s." % (fdir+fName+'.hdf5'))
            
    # Tips to open and read the hdf5 file in python:
    #    g=h5py.File('filename.hdf5','r')
    # For arrays:
    #    data=g["data"][:]      The whole array is returned in variable "data"
    # For dictionaries containing arrays:
    #    t=g["times"]
    #    t["abovecloud"][:]
    # For dictionaries containing strings:
    #    desc=g["desc"]
    #    desc2=dict()
    #    for key in desc.attrs:
    #        desc2[key]=desc.attrs[key]
    # For lists, now np.strings:
    #    dttl=g["dttl"]
    #    L=str(dttl[...]).replace(" u'","").replace(" '","").replace("'","").strip("[").strip("]").split(",")


