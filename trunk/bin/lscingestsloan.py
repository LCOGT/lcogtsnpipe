#!/science/supernova/anaconda/bin/python

if __name__ == '__main__':
    import lsc
    from lsc import conn
    from lsc import readkey3, readhdr
    import argparse

    parser = argparse.ArgumentParser(description='Download and ingest SDSS or PS1 templates')
    parser.add_argument('images', nargs='+', help='LCOGT images to use for field of view and filter')
    parser.add_argument('--type', default='sloan', help='Survey to download from', choices=['sloan', 'ps1'])
    parser.add_argument('--ps1frames', default='')
    args = parser.parse_args()
    imglist = args.images
    imgtype = args.type
    ps1frames = args.ps1frames

    if imgtype =='sloan':
        for img in imglist:
            image0, varimg = lsc.sloanimage(img,'sloan','')
    elif imgtype =='ps1':
        print "WARNING: PS1 ingestion works at the moment with single object and filter\n "
        print "please, do not provide multiple objects and filter in the same query"
#        if not ps1frames:
#            sys.exit('ERROR: you need to provide the PS1 files')
#        else:
        if ps1frames:
            frames = np.genfromtxt(ps1frames,str)
        else:
            frames=''
        for img in imglist:
            image0, varimg = lsc.sloanimage(img,'ps1',frames)
    else:
        image0=''
        print 'add here ingestion of different images (DES)'

    if image0:
        lsc.mysqldef.ingestdata('tar', '', ['20121212'], False, 'oracproc', '', [image0])
        lsc.mysqldef.ingestredu([image0], True, 'photlco')
        _path = lsc.mysqldef.query(['select filepath from photlco where filename = "'+str(image0)+'"'],conn)
        if _path and varimg:
            print _path[0]['filepath']
            print varimg
            os.system('cp '+varimg+ ' '+_path[0]['filepath'])
            print 'cp '+varimg+ ' '+_path[0]['filepath']
        lsc.mysqldef.updatevalue('photlco','filetype',4,image0,'lcogt2','filename')        
