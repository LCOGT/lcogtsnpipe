
#######################################################################
def dbConnect(lhost, luser, lpasswd, ldb):
   import sys
   import MySQLdb,os,string
   try:
      conn = MySQLdb.connect (host = lhost,
                              user = luser,
                            passwd = lpasswd,
                                db = ldb)
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      sys.exit (1)
   return conn

def getconnection(site):
   import lsc
   connection = {'lcogt2':{}, 'iptf':{}}
   connection['lcogt2']['database'] = lsc.util.readpass['database'] 
   connection['lcogt2']['hostname'] = lsc.util.readpass['hostname'] 
   connection['lcogt2']['username'] = lsc.util.readpass['mysqluser'] 
   connection['lcogt2']['passwd'] = lsc.util.readpass['mysqlpasswd'] 

   connection['iptf']['database'] = lsc.util.readpass['ptfdatabase'] 
   connection['iptf']['hostname'] = lsc.util.readpass['ptfhost'] 
   connection['iptf']['username'] = lsc.util.readpass['ptfuser'] 
   connection['iptf']['passwd'] = lsc.util.readpass['ptfpasswd'] 

   return  connection[site]['hostname'],connection[site]['username'],connection[site]['passwd'],connection[site]['database']

########################################################################

def getmissing(conn, epoch0, epoch2,telescope,datatable='photlco'):
   import sys
   import lsc
   import MySQLdb,os,string
   print epoch0, epoch2,telescope
   try:
      cursor = conn.cursor (MySQLdb.cursors.DictCursor)
      if telescope =='all':
         if epoch2:
            print "select raw.filename, raw.objname from photlcoraw as raw where "+\
                " raw.dayobs < "+str(epoch2)+" and raw.dayobs >= "+str(epoch0)+\
                " and NOT EXISTS(select * from "+str(datatable)+" as redu where raw.filename = redu.filename)"
            cursor.execute ("select raw.filename, raw.objname from photlcoraw as raw where "+\
                               " raw.dayobs < "+str(epoch2)+" and raw.dayobs >= "+str(epoch0)+\
                               " and NOT EXISTS(select * from "+str(datatable)+" as redu where raw.filename = redu.filename)")
         else:
            cursor.execute ("select raw.filename, raw.objname from photlcoraw raw where "+\
                               " raw.dayobs = "+str(epoch0)+\
                               " and NOT EXISTS(select * from "+str(datatable)+" redu where raw.filename = redu.filename)")
      else:
         fntel = telescope.replace('-', '') # 1m0-01 (input) --> 1m001 (in filename)
         if epoch2:  
            print "select raw.filename, raw.objname from photlcoraw raw where (raw.filename like '%"+fntel+"%'"+\
                               " or raw.telescope = '"+telescope+"') and raw.dayobs < "+str(epoch2)+" and raw.dayobs >= "+str(epoch0)+\
                               " and NOT EXISTS(select * from "+str(datatable)+" redu where raw.filename = redu.filename)"
            cursor.execute ("select raw.filename, raw.objname from photlcoraw raw where (raw.filename like '%"+fntel+"%'"+\
                               " or raw.telescope = '"+telescope+"') and raw.dayobs < "+str(epoch2)+" and raw.dayobs >= "+str(epoch0)+\
                               " and NOT EXISTS(select * from "+str(datatable)+" redu where raw.filename = redu.filename)")
         else:
            print "select raw.filename, raw.objname from photlcoraw raw where (raw.filename like '%"+fntel+"%'"+\
                               " or raw.telescope = '"+telescope+"') and raw.dayobs = "+str(epoch0)+\
                               " and NOT EXISTS(select * from "+str(datatable)+" redu where raw.filename = redu.filename)"
            cursor.execute ("select raw.filename, raw.objname from photlcoraw raw where (raw.filename like '%"+fntel+"%'"+\
                               " or raw.telescope = '"+telescope+"') and raw.dayobs = "+str(epoch0)+\
                               " and NOT EXISTS(select * from "+str(datatable)+" redu where raw.filename = redu.filename)")
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      sys.exit (1)
   return resultSet

def getfromdataraw(conn, table, column, value,column2='*'):
   import sys
   import MySQLdb,os,string
   try:
      cursor = conn.cursor (MySQLdb.cursors.DictCursor)
      cursor.execute ("select "+column2+" from "+str(table)+" where "+column+"="+"'"+value+"'")
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      sys.exit (1)
   return resultSet


def getlistfromraw(conn, table, column, value1, value2, column2='*', telescope='all'):
    import MySQLdb
    cursor = conn.cursor(MySQLdb.cursors.DictCursor)
    if not value2:
        value2 = value1
    if telescope=='all':
        cursor.execute("select "+column2+" from "+str(table)+" where "+column+"<="+"'"+value2+"' and "+column+">="+"'"+value1+"'")
    else:
        fntel = telescope.replace('-', '')  # 1m0-01 (input) --> 1m001 (in filename)
        cursor.execute("select "+column2+" from "+str(table)+" where "+column+"<="+"'"+value2+"' and "+column+">="+"'"+value1+"' and (filename like '%"+fntel+"%' or telescope='"+telescope+"')")
    resultSet = cursor.fetchall()
    cursor.close()
    return resultSet

###########################################################################

def updatevalue(table,column,value,filename,connection='lcogt2',filename0='filename'):
   import sys
   import MySQLdb,os,string
   import lsc

   hostname, username, passwd, database=lsc.mysqldef.getconnection(connection)
   conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
   try:
      cursor = conn.cursor (MySQLdb.cursors.DictCursor)
      if isinstance(column, list):
         columns = column
      else:
         columns = [column]
      if isinstance(value, list):
         values = ["'{}'".format(val) if isinstance(val, str) else val for val in value]
      elif isinstance(value, str):
         values = ["'{}'".format(value)]
      else:
         values = [value]
      column_equals_value = ', '.join(['{}={}'.format(col, val) for col, val in zip(columns, values)])
      query = "UPDATE {} SET {} WHERE {}='{}'".format(table, column_equals_value, filename0, filename)
      print query
      cursor.execute(query)
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      conn.commit()
      cursor.close ()
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])

###########################################################################

def insert_values(conn,table,values):
    import sys,string,os,re,MySQLdb,os,string,datetime

    datecreated_tables = ['atels','groups','iaunames','instruments','notes',
                          'obsrequests','photlco','photlcoraw','photpairing',
                          'spec','speclcoraw','targetnames','targets',
                          'telescopes','useractionlog','users']
    if 'datecreated' not in values and table in datecreated_tables:
        values['datecreated'] = str(datetime.datetime.utcnow())

    def dictValuePad(key):
        return '%(' + str(key) + ')s'

    def insertFromDict(table, dicto):
        """Take dictionary object dict and produce sql for 
        inserting it into the named table"""
        cleandict = {key: val for key, val in dicto.items() if val not in ['NaN', 'UNKNOWN', 'N/A', None, '']}

        sql = 'INSERT INTO ' + table
        sql += ' ('
        sql += ', '.join(cleandict)
        sql += ') VALUES ('
        sql += ', '.join(map(dictValuePad, cleandict))
        sql += ');'
        return sql

    sql = insertFromDict(table, values)
    try:
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        cursor.execute(sql, values)
        resultSet = cursor.fetchall ()
        if cursor.rowcount == 0:
            pass
        conn.commit()
        cursor.close ()
    except MySQLdb.Error, e:
        print "Error %d: %s" % (e.args[0], e.args[1])

########################################################################

def guess_instrument_type(name):
    prefix = name[:2]
    if prefix == 'fl' or prefix == 'fa':
        insttype = 'Sinistro'
    elif prefix == 'kb':
        insttype = 'SBIG'
    elif prefix == 'fs':
        insttype = 'Spectral'
    elif prefix == 'ep':
        insttype = 'MUSCAT'
    else:
        insttype = None
    return insttype

def ingestredu(imglist,force='no',dataredutable='photlco',filetype=1):
   import string,re,os,sys
   import lsc
   from lsc.util import readkey3, readhdr
   from datetime import datetime

   hostname, username, passwd, database=lsc.mysqldef.getconnection('lcogt2')
   conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)

   for fullpath in imglist:
      print(fullpath)
      path, img = os.path.split(fullpath)
      path += '/'

######  for external users ingesting directly in lcophot ######      
      if img[-3:] == '.fz':
         os.system('funpack -D ' + fullpath)
         img = img[:-3]
         fullpath = fullpath[:-3]
         if os.path.isfile(fullpath + '.fz'):
            os.remove(fullpath + '.fz')
############################################################

      exist=lsc.mysqldef.getfromdataraw(conn,dataredutable,'filename', string.split(img,'/')[-1],column2='filename')
      exist2=lsc.mysqldef.getfromdataraw(conn,'photlcoraw','filename', string.split(img,'/')[-1],column2='filename, groupidcode')
      if exist2:
         print exist2
         _groupidcode=exist2[0]['groupidcode']
      else:
         _groupidcode=''

      if exist:
         if force=='yes':
            print img,database
            lsc.mysqldef.deleteredufromarchive(string.split(img,'/')[-1],dataredutable)
            print 'delete line from '+str(database)
            exist = False
            
      if not exist or force =='update':
         hdr = readhdr(fullpath)
         _targetid = lsc.mysqldef.targimg(fullpath)
         try:
            _tracknumber=int(readkey3(hdr,'TRACKNUM'))
         except:
            _tracknumber=0
         if hdr.get('TELESCOP'):  
            _tel=hdr.get('TELESCOP')
         else: 
            _tel=''
         if _tel in ['Faulkes Telescope South','fts']:  
            _tel='2m0-02'
         elif _tel in ['Faulkes Telescope North','ftn']: 
            _tel='2m0-01'
         _inst=hdr.get('instrume')
         dictionary={'dateobs':readkey3(hdr,'date-obs'),'dayobs':readkey3(hdr,'DAY-OBS'),'filename':img,'filepath':path,'filetype': int(filetype), 'targetid':_targetid,\
                     'exptime':readkey3(hdr,'exptime'), 'filter':readkey3(hdr,'filter'),'mjd':readkey3(hdr,'mjd'),'tracknumber':_tracknumber,'groupidcode':_groupidcode,\
                     'telescope':_tel,'airmass':readkey3(hdr,'airmass'),'objname':readkey3(hdr,'object'),'ut':readkey3(hdr,'ut'),\
                     'wcs':readkey3(hdr,'wcserr'),'instrument':_inst,'ra0':readkey3(hdr,'RA'),'dec0':readkey3(hdr,'DEC')}

         _telid=lsc.mysqldef.getfromdataraw(conn,'telescopes','name',_tel,column2='id')
         if not _telid:
           print 'Telescope ',_tel,' not recognized.  Adding to telescopes table.'
           # the short name is needed to calibrate the magnitude with extinction
           lsc.mysqldef.insert_values(conn, 'telescopes', {'name': _tel, 'shortname': readkey3(hdr,'SITEID')})
           
           _telid=lsc.mysqldef.getfromdataraw(conn,'telescopes','name',_tel,column2='id')
         telid=_telid[0]['id']
         dictionary['telescopeid']=str(telid)

         _instid=lsc.mysqldef.getfromdataraw(conn,'instruments','name',_inst,column2='id')
         if not _instid:
           print 'Instrument ',_inst,' not recognized.  Adding to instruments table.'
           lsc.mysqldef.insert_values(conn, 'instruments', {'name': _inst, 'type': guess_instrument_type(_inst)})
           _instid=lsc.mysqldef.getfromdataraw(conn,'instruments','name',_inst,column2='id')
         instid=_instid[0]['id']
         dictionary['instrumentid']=str(instid)
         dictionary['lastunpacked'] = str(datetime.utcnow())

         print dictionary
         print 'insert reduced'

         # we need to update the connection before quering again the database
         hostname, username, passwd, database=lsc.mysqldef.getconnection('lcogt2')
         conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
         
         ggg=lsc.mysqldef.getfromdataraw(conn, dataredutable, 'filename',str(img), '*')
         if not ggg:
            lsc.mysqldef.insert_values(conn,dataredutable,dictionary)
         else:
            for voce in dictionary:
               lsc.mysqldef.updatevalue(dataredutable,voce,dictionary[voce],string.split(img,'/')[-1])
      else:
         print 'already ingested'

###############################################################################################################################

def getvaluefromarchive(table,column,value,column2):
   import sys
   import MySQLdb,os,string
   import lsc
   hostname, username, passwd, database=lsc.mysqldef.getconnection('lcogt2')
   conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
   resultSet=lsc.mysqldef.getfromdataraw(conn, table, column, value,column2)
   if resultSet:
      result=resultSet
   else:
      result=[]
   return result


def deleteredufromarchive(filename,archive='photlco',column='filename'):
   import sys
   import MySQLdb,os,string
   import lsc
   hostname, username, passwd, database=lsc.mysqldef.getconnection('lcogt2')
   conn = lsc.mysqldef.dbConnect(hostname, username, passwd, database)
#######
   try:
      cursor = conn.cursor (MySQLdb.cursors.DictCursor)
      cursor.execute ("delete  from "+str(archive)+" where "+str(column)+"="+"'"+filename+"'")
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      conn.commit()
      cursor.close ()
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      sys.exit (1)
   return resultSet

###################################################################################################################
def updateDatabase(tarfile):
   import string,os,re,math,sys,shutil,glob,socket,pickle
   import MySQLdb
   host=socket.gethostname()
   import ntt
   from ntt.util import readkey3,readhdr
   from mysqldef import updatevalue, getvaluefromarchive, dbConnect,getfromdataraw,ingestredu
   os.system('tar -zxvf '+str(tarfile))
   pklfile=re.sub('tar.gz','pkl',tarfile)
   if glob.glob(pklfile):
      pkl_file = open(pklfile, 'rb')
      mydict2 = pickle.load(pkl_file)
      pkl_file.close()
      fitsfile=[]
      for i in mydict2:
         try:
            command=mydict2[i]['command']
            if command=='ingestredu':
               lista=mydict2[i]['lista']
               instrument=mydict2[i]['instument']
               output=mydict2[i]['output']
               ingestredu(lista)
               if output not in fitsfile: fitsfile.append(output)
            elif command=='getvaluefromarchive':
               table=mydict2[i]['table']
               column=mydict2[i]['column']
               value=mydict2[i]['value']
               _input=mydict2[i]['input']
               hh=getvaluefromarchive(table,column,_input,value)
            elif command=='updatevalue':
               table=mydict2[i]['table']
               column=mydict2[i]['column']
               value=mydict2[i]['value']
               _input=mydict2[i]['input']
               if value!='voce':
                  updatevalue(table,column,value,_input)
               else:
                  if column=='voce':
                     for voce in hh:
                        if voce!='id' and voce!='filename':
                           updatevalue(table,voce,hh[voce],_input)
                  else:
                     updatevalue(table,column,hh[column],_input)
            else:
               print 'warning: command  not recognise'
         except:
            print '#####################'
            print mydict2[i]
            print 'problems'
      print fitsfile
      dir1='/data/obsdata/y20'+string.split(tarfile,'_')[2][-2:]+'/'
      dir2='/data/obsdata/y20'+string.split(tarfile,'_')[2][-2:]+'/'+string.split(tarfile,'_')[2]+'/'
      if not os.path.isdir(dir1):  os.mkdir(dir1)
      if not os.path.isdir(dir2):  os.mkdir(dir2)
      for img in fitsfile:
         arcfile=readkey3(readhdr(img),'ARCFILE')
         print arcfile
         try:
            bbbb=getvaluefromarchive('datarawNTT','filename',arcfile,'*')
         except:  bbbb=''
         if bbbb:
            directoryraw=bbbb['filepath']
            directoryreduced=re.sub('raw/','',directoryraw)+'reduced/'
         else:     
            dir3=dir2+string.split(tarfile,'_')[3]+'_'+string.split(tarfile,'_')[4]+'/'
            if not os.path.isdir(dir3):  os.mkdir(dir3)
            directoryreduced=dir3+'reduced/'

         if not os.path.isdir(directoryreduced):        os.mkdir(directoryreduced)
         os.system('cp '+img+' '+directoryreduced)
         updatevalue('redulogNTT','filepath',directoryreduced,img)
         updatevalue('datareduNTT','filepath',directoryreduced,img)

##############################################################################

def getfromcoordinate(conn, table, ra0, dec0,distance):
   import sys
   import MySQLdb,os,string
   #  this is acually not needed
   if table=='targets':
      ra1='ra0'
      dec1='dec0'
   elif table=='photlcoraw':
      ra1='ra0'
      dec1='dec0'
   ####
   try:
      cursor = conn.cursor (MySQLdb.cursors.DictCursor)
      command=["set @sc = pi()/180","set @ra = "+str(ra0), "set @dec = "+str(dec0),"set @distance = "+str(distance),"SELECT *,abs(2*asin( sqrt( sin((a.dec0-@dec)*@sc/2)*sin((a.dec0-@dec)*@sc/2) + cos(a.dec0*@sc)*cos(@dec*@sc)*sin((a.ra0-@ra)*@sc/2)*sin((a.ra0-@ra)*@sc/2.0) )))*180/pi() as hsine FROM "+str(table)+" as a HAVING hsine<@distance order by a.ra0 desc"]
      for ccc in command:
         cursor.execute (ccc)
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      sys.exit (1)
   return resultSet

###################################################################################

def targimg(img='', hdrt=None):
    import lsc
    from lsc.util import readkey3,readhdr
    from lsc.mysqldef import getfromcoordinate
    from lsc import conn
    import string
    import os, json, requests
    _targetid=''
    _group=''
    if hdrt is None:
        hdrt=lsc.util.readhdr(img)
    
    _ra=lsc.util.readkey3(hdrt,'CAT-RA')
    _dec=lsc.util.readkey3(hdrt,'CAT-DEC')

    if _ra is None or _dec is None:
        # No CAT RA or dec, so send a warning message to slack and return exception
        # Send Slack message
        post_url = os.environ['SLACK_CHANNEL_WEBHOOK']
        payload = {'text': 'CAT-RA and CAT-DEC could not be found for {}'.format(img)}
        json_data = json.dumps(payload)
        headers = {'Content-Type': 'application/json'}
        response = requests.post(post_url, data=json_data.encode('ascii'), headers=headers)

        # Raise exception so pipeline moves on to ingesting the next image
        raise Exception ('No CAT-RA or CAT-DEC could be found for {}'.format(img))

    _object=lsc.util.readkey3(hdrt,'object')
    if ':' in str(_ra):        
       _ra,_dec=lsc.deg2HMS(_ra,_dec)

    #############  define groupid ################
    aaa=lsc.mysqldef.query(['select idcode, groupidcode from programs'],conn)
    grname=[i['idcode'] for i in aaa]
    gr=[i['groupidcode'] for i in aaa]
    if lsc.util.readkey3(hdrt,'propid') in grname:
       if gr[grname.index(lsc.util.readkey3(hdrt,'propid'))] not in [1,None]:
          _group = gr[grname.index(lsc.util.readkey3(hdrt,'propid'))]
       else:
          _group=32769
    else:
          _group=32769
    ##############################################

    ########  define targetid  ##################
    _targetid=lsc.mysqldef.gettargetid(_object,'','',conn,.01,False)
    if not _targetid:
        print '# no target with this name '+_object
        _targetid=lsc.mysqldef.gettargetid('',_ra,_dec,conn,.01,False)
        if _targetid:
            print '# target at this coordinate with a different name, add name '+str(_ra)+' '+str(_dec)
            dictionary1={'name':_object,'targetid':_targetid,'groupidcode':_group}
            lsc.mysqldef.insert_values(conn,'targetnames',dictionary1)

            bb2=lsc.mysqldef.query(['select id from targetnames where name = "'+_object+'"'],conn)
            dictionary3 = {'userid':67, 'tablemodified': 'targetnames', 'idmodified': bb2[0]['id'],
                          'columnmodified': 'New Row', 'newvalue': 'Multiple'}
            lsc.mysqldef.insert_values(conn,'useractionlog',dictionary3)
        else:
            print 'not found targetid with ra and dec '+str(_ra)+' '+str(_dec)

    else:
        print 'found name= '+_object+'  targetid= '+str(_targetid)
    ##############################################

    if not _targetid:
          # no target
          print 'add new target '+str(_ra)+' '+str(_dec)+' '+_object
          dictionary={'ra0':_ra,'dec0':_dec}
          lsc.mysqldef.insert_values(conn,'targets',dictionary)
          bb=lsc.mysqldef.getfromcoordinate(conn, 'targets', _ra, _dec,.000156)

          dictionary2 = {'userid':67, 'tablemodified': 'targets', 'idmodified': bb[0]['id'],
                         'columnmodified': 'New Row', 'newvalue': 'Multiple'}
          lsc.mysqldef.insert_values(conn,'useractionlog',dictionary2)

          dictionary1 = {'name':_object,'targetid':bb[0]['id'],'groupidcode':_group}
          lsc.mysqldef.insert_values(conn,'targetnames',dictionary1)

          bb2=lsc.mysqldef.query(['select id from targetnames where name = "'+_object+'"'],conn)
          dictionary3 = {'userid':67, 'tablemodified': 'targetnames', 'idmodified': bb2[0]['id'],
                         'columnmodified': 'New Row', 'newvalue': 'Multiple'}
          lsc.mysqldef.insert_values(conn,'useractionlog',dictionary3)


          _targetid=bb[0]['id']
          print '\n add a target '+str(_ra)+' '+str(_dec)+' '+str(bb[0]['id'])

    if _targetid and _group:
       cc=lsc.mysqldef.getfromdataraw(conn,'permissionlog','targetid', str(_targetid),column2='groupname')
       if len(cc)==0:
          _JDn=lsc.mysqldef.JDnow()
          print img
          dictionary2={'targetid':_targetid,'jd':_JDn,'groupname':_group}
          lsc.mysqldef.insert_values(conn,'permissionlog',dictionary2)
    return _targetid

#################################################################################
def getsky(data):
  """
  Determine the sky parameters for a FITS data extension.

  data -- array holding the image data
  """
  from numpy import random

  # maximum number of interations for mean,std loop
  maxiter = 30

  # maximum number of data points to sample
  maxsample = 10000

  # size of the array
  ny,nx = data.shape

  # how many sampels should we take?
  if data.size > maxsample:
    nsample = maxsample
  else:
    nsample = data.size

  # create sample indicies
  xs = random.uniform(low=0, high=nx, size=nsample).astype('L')
  ys = random.uniform(low=0, high=ny, size=nsample).astype('L')

  # sample the data
  sample = data[ys,xs].copy()
  sample = sample.reshape(nsample)

  # determine the clipped mean and standard deviation
  mean = sample.mean()
  std = sample.std()
  oldsize = 0
  niter = 0
  while oldsize != sample.size and niter < maxiter:
    niter += 1
    oldsize = sample.size
    wok = (sample < mean + 3*std)
    sample = sample[wok]
    wok = (sample > mean - 3*std)
    sample = sample[wok]
    mean = sample.mean()
    std = sample.std()
 
  return mean,std

################################################################################

def getlike(conn, table, column, value,column2='*'):
   import sys
   import MySQLdb,os,string
   try:
      cursor = conn.cursor (MySQLdb.cursors.DictCursor)
      cursor.execute ("select "+column2+" from "+str(table)+" where "+column+" like "+"'%"+value+"%'")
      resultSet = cursor.fetchall ()
      if cursor.rowcount == 0:
         pass
      cursor.close ()
   except MySQLdb.Error, e:
      print "Error %d: %s" % (e.args[0], e.args[1])
      sys.exit (1)
   return resultSet

##################################################################

def query(command,conn):
   import MySQLdb,os,string
   lista=''
   #from lsc import conn
   try:
        cursor = conn.cursor (MySQLdb.cursors.DictCursor)
        for i in command:
            cursor.execute (i)
            lista = cursor.fetchall ()
            if cursor.rowcount == 0:
                pass
        conn.commit()
        cursor.close ()
   except MySQLdb.Error, e: 
        print "Error %d: %s" % (e.args[0], e.args[1])
   return lista

###############################
def JDnow(datenow='',verbose=False):
   import datetime
   import time
   _JD0=2455927.5
   if not datenow:
      datenow=datetime.datetime(time.gmtime().tm_year, time.gmtime().tm_mon, time.gmtime().tm_mday, time.gmtime().tm_hour, time.gmtime().tm_min, time.gmtime().tm_sec)
   _JDtoday=_JD0+(datenow-datetime.datetime(2012, 01, 01,00,00,00)).seconds/(3600.*24)+\
             (datenow-datetime.datetime(2012, 01, 01,00,00,00)).days
   if verbose: print 'JD= '+str(_JDtoday)
   return _JDtoday
###################################

###############################
def MJDnow(datenow='',verbose=False):
   import datetime
   import time
   _JD0=55927.
   if not datenow:
      datenow=datetime.datetime(time.gmtime().tm_year, time.gmtime().tm_mon, time.gmtime().tm_mday, time.gmtime().tm_hour, time.gmtime().tm_min, time.gmtime().tm_sec)
   _JDtoday=_JD0+(datenow-datetime.datetime(2012, 01, 01,00,00,00)).seconds/(3600.*24)+\
             (datenow-datetime.datetime(2012, 01, 01,00,00,00)).days
   if verbose: print 'JD= '+str(_JDtoday)
   return _JDtoday
###################################

def gettargetid(_name,_ra,_dec,conn,_radius=.01,verbose=False):
   import lsc
   from numpy import argmin
   if _name:
        command=['select distinct(r.name), r.targetid, t.id, t.ra0, t.dec0 from targetnames as r join targets as t where r.name like "%'+\
                 _name.replace(' ', '%') +'" and t.id=r.targetid']
        lista=lsc.mysqldef.query(command,conn)
   elif _ra and _dec:
      if ':' in  str(_ra):
         _ra,_dec=lsc.deg2HMS(_ra,_dec) 
      lista=lsc.mysqldef.getfromcoordinate(conn, 'targets', _ra, _dec,_radius)
      if verbose:
         print _ra,_dec,_radius
         print lista
   if lista:
        ll0={}
        for jj in lista[0].keys(): ll0[jj]=[]
        for i in range(0,len(lista)):
            for jj in lista[0].keys(): ll0[jj].append(lista[i][jj])
            ###########################################        
        if len(set(ll0['id']))==1:
            _targetid=ll0['id'][0]
        elif  len(set(ll0['id']))>1:
           if 'hsine' in ll0.keys():
              _targetid=ll0['id'][argmin(ll0['hsine'])]
           else:
              _targetid=''
        else:
            _targetid=''
   else:
        _targetid=''
        if verbose:   
            print 'no objects'
   return _targetid

###############################################################
def get_snex_uid(interactive=True, return_fullname=False):
    from getpass import getuser
    from lsc import conn
    usersdict = query(['select id, firstname, lastname from users where name="{}"'.format(getuser())], conn)
    if usersdict: # if your UNIX username matches your SNEx username
        snex_uid = usersdict[0]['id']
        fullname = usersdict[0]['firstname'] + ' ' + usersdict[0]['lastname']
    else:
        print 'Your username is not associated with a SNEx account.'
        snex_uid = None
        fullname = None
        if interactive:
            snex_user = raw_input('If you have a SNEx username, input it here. Otherwise, press enter. ')
            if snex_user:
                usersdict = query(['select id, firstname, lastname from users where name="' + snex_user + '"'], conn)
                if usersdict:
                  snex_uid = usersdict[0]['id']
                  fullname = usersdict[0]['firstname'] + ' ' + usersdict[0]['lastname']
                else:
                  print 'SNEx username not found.'
    if return_fullname: return snex_uid, fullname
    else: return snex_uid
