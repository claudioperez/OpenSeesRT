
add_library(SuperLU_MT)

target_compile_definitions(SuperLU_MT PRIVATE _PTHREAD)

target_sources(SuperLU_MT
  PUBLIC
#   util.h
#   Cnames.h
#   supermatrix.h
#   colamd.h
#   machines.h
#   pdsp_defs.h
#   pxgstrf_synch.h

  PRIVATE
  # ALLAUX
    superlu_timer.c dclock.c sp_ienv.c lsame.c xerbla.c
    util.c pmemory.c qrnzcnt.c dmyblas2.c await.c
    get_perm_c.c mmd.c colamd.c sp_coletree.c sp_colorder.c
    pxgstrf_relax_snode.c
    pxgstrf_synch.c
    pxgstrf_scheduler.c
    pxgstrf_mark_busy_descends.c
    pxgstrf_pruneL.c
    pxgstrf_super_bnd_dfs.c

  # DZAUX
    dlamch.c

  # DLUSRC
    dreadhb.c dmatgen.c
    pdgssv.c pdgssvx.c dgstrs.c dgsrfs.c
    dgscon.c dlacon.c dlangs.c dgsequ.c dlaqgs.c dpivotgrowth.c
    pdmemory.c pdutil.c
    pdgstrf.c pdgstrf_init.c pdgstrf_finalize.c
    pdgstrf_thread.c pdgstrf_thread_init.c pdgstrf_thread_finalize.c
    pdgstrf_factor_snode.c pdgstrf_snode_dfs.c pdgstrf_snode_bmod.c
    pdgstrf_panel_dfs.c pdgstrf_panel_bmod.c
    pdgstrf_copy_to_ucol.c pdgstrf_pivotL.c
    pdgstrf_column_dfs.c pdgstrf_column_bmod.c
    pdgstrf_bmod1D.c pdgstrf_bmod2D.c
    pdgstrf_bmod1D_mv2.c pdgstrf_bmod2D_mv2.c
    dsp_blas2.c dsp_blas3.c
)
