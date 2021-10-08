#!/usr/bin/env python3
"""
    Imports SAdLSA data into an sqlite3 table
"""



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='SAdLSA data importer')
    parser.add_argument('--sadlsa-dir', required=True,
                        help='directory containing SAdLSA alignment and score filed got a single protein')
    parser.add_argument('--database', required=True,
                        help='sqlite3 in which to insert EC data')

    args = parser.parse_args()

    sadlsa_dir = Path(args.sadlsa_dir)
    if not sadlsa_dir.exists():
        logging.critical(f'{args.sadlsa_dir} does not exist ... exiting')
        sys.exit(1)

    logging.info(f'Ingesting {args.sadlsa_dir}')

    ec_df = score_to_sqlite3(sadlsa_dir)

    write_df_to_db(ec_df, args.database)

    logging.info(f'Done.  Added {len(ec_df)} entries to'
                 f' {args.database} table {EC_TABLE}.')
