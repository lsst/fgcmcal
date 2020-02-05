from lsst.pipe.tasks.colorterms import Colorterm, ColortermDict

config.data = {
    "sdss*": ColortermDict(data={
            'g': Colorterm(primary='g', secondary='r',
                           c0=-0.00816446, c1=-0.08366937, c2=-0.00726883),
            'r': Colorterm(primary='r', secondary='i',
                           c0=0.0013181, c1=0.01284177, c2=-0.03068248),
            'i': Colorterm(primary='i', secondary='z',
                           c0=0.00130204, c1=-0.16922042, c2=-0.01374245)
    }),
}
