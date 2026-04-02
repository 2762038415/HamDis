from hammeth.scripts.hamdist_core import run_hamdist as hamdist_core


def run_hamdist(args):
    hamdist_core(
        pat_dir=args.i,
        output_dir=args.o,
        batch_size=args.bs,
        batch_num=args.bn,
        m_process=args.m,
        bin_num=args.n,
        region_dir=args.r,
        cpg_file_dir=args.cp,
        genome=args.g,
        chr_list=args.cl,
        binsize=getattr(args, "binsize", 5000),
        seed=getattr(args, "seed", None),
        min_hm_count=getattr(args, "min_hm_count", 5),
        min_cpg=getattr(args, "min_cpg", 2),
        max_dis=getattr(args, "max_dis", 600),
        max_possible_read=getattr(args, "max_possible_read", 3000),
        basedon0=getattr(args, "basedon0", 0),
        vmr=getattr(args, "vmr", "hmtop"),
        vmr_bed=getattr(args, "vmr_bed", None),
    )
