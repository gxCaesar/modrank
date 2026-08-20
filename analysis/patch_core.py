import pathlib
p = pathlib.Path("utils/core_utils.py")
s = p.read_text()
old = "        model_dict = {'omic_sizes': args.omic_sizes, 'num_classes': args.n_classes}\n\n        if args.use_nystrom:"
new = (
    "        # UPSTREAM BUG, wired through rather than worked around. SurvPath.__init__ accepts\n"
    "        # wsi_embedding_dim and defaults it to 1024, but this model_dict never passes it, so\n"
    "        # --encoding_dim is silently ignored for the survpath modality and reaches only the\n"
    "        # results-directory name (general_utils.py:135). With the 768-d CTransPath features the\n"
    "        # SurvPath paper itself specifies, the run dies with\n"
    "        #   RuntimeError: mat1 and mat2 shapes cannot be multiplied (4096x768 and 1024x256)\n"
    "        # Passing the existing flag through is what makes the released code do what its own\n"
    "        # README describes; it does not alter the method.\n"
    "        model_dict = {'omic_sizes': args.omic_sizes, 'num_classes': args.n_classes,\n"
    "                      'wsi_embedding_dim': args.encoding_dim}\n\n        if args.use_nystrom:"
)
assert s.count(old) == 1, "anchor count = %d" % s.count(old)
p.write_text(s.replace(old, new))
print("patched utils/core_utils.py")
