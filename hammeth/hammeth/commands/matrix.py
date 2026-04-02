from hammeth.scripts.matrix import run_matrix_pipeline as matrix_core


def run_matrix(args):
    matrix_core(
        input_dir=args.i,
        region_bed=args.region_bed,
        tmp_dir=args.tmp_dir,
        out_long=args.out_long,
        out_matrix=args.out_matrix,
        memmap_file=args.memmap_file,
        nproc=args.nproc,
        batch_size=args.batch_size,
        overwrite=getattr(args, "overwrite", False),
        keep_tmp=getattr(args, "keep_tmp", False),
    )
