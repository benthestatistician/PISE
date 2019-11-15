Setup for this, using 
    R 3.4.1 / x86_64-apple-darwin15.6.0 (64-bit)

and also:
- R 3.0.2/ x86_64-pc-linux-gnu (64-bit)
- R 3.6.1/  x86_64-apple-darwin15.6.0 (64-bit)

```r
devtools::install_github('ramnathv/slidify@1dd41a3')
devtools::install_github('ramnathv/slidifyLibraries@dbd065f')
```

Then added directories and file
      assets/layouts/twocolleftwider.html
      

And then I was able to generate from within this directory, using 

    slidify::slidify("index.Rmd")
