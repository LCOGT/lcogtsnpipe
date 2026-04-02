#!/usr/bin/env python
import lsc
import os
from LCOGTingest import db_ingest
from argparse import ArgumentParser
import logging

logger = logging.getLogger()


if __name__ == "__main__":
    parser = ArgumentParser(description='Ingest data from the Supernova Exchange')
    parser.add_argument("frames", nargs='+', help="images to ingest")
    parser.add_argument("-F", "--force", action="store_true", help="reingest files even if they already exist")
    args = parser.parse_args()

    logger.info(f'Ingesting {len(args.frames)} frames')

    fullpaths = [os.path.abspath(frame) for frame in args.frames]
    for frame in fullpaths:
        dbdict = db_ingest(*os.path.split(frame), force=args.force)
    lsc.mysqldef.ingestredu(fullpaths, force='yes' if args.force else 'no')  # ingest new data into photlco
