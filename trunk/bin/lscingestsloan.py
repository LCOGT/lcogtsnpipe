#!/usr/bin/env python

if __name__ == '__main__':
    import lsc
    from lsc import conn
    from lsc import readkey3, readhdr
    import argparse
    import os
    from astropy.io import fits

    parser = argparse.ArgumentParser(description='Download and ingest SDSS or PS1 templates')
    parser.add_argument('images', nargs='+', help='LCOGT images to use for field of view and filter')
    parser.add_argument('--type', default='sloan', help='Survey to download from', choices=['sloan', 'ps1'])
    parser.add_argument('--ps1frames', default='')
    parser.add_argument('--show', action='store_true', help='show the assembled, swarped SDSS template in DS9')
    parser.add_argument('-F', '--force', action='store_true', help='redownload SDSS images if already present in CWD')
    args = parser.parse_args()
    imglist = args.images
    imgtype = args.type
    ps1frames = args.ps1frames

    if imgtype =='sloan':
        for img in imglist:
            image0, varimg = lsc.sloanimage(img,'sloan','', args.show, args.force)
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
            image0, varimg = lsc.sloanimage(img,'ps1',frames, args.show)
    else:
        image0=''
        print 'add here ingestion of different images (DES)'

    if image0:
        hdr = fits.getheader(image0)
        filepath = lsc.util.workdirectory + 'data/extdata/' + hdr['DAYOBS'] + '/'
        os.system('mv -v {} {} {}'.format(image0, varimg, filepath))
        lsc.LCOGTingest.db_ingest(filepath, image0, force=args.force)
        lsc.mysqldef.ingestredu([filepath + image0], 'yes' if args.force else 'no', 'photlco', filetype=4)
