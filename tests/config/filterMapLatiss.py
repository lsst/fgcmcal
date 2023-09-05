for source, target in [
        ("SDSSg_65mm~empty", "g"),
        ("SDSSr_65mm~empty", "r"),
        ("SDSSi_65mm~empty", "i"),
]:
    config.filterMap[source] = target
