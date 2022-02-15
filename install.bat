echo y | pyinstaller --onedir --icon="pi.ico" PyCalc.pyw
xcopy pi.ico dist\PyCalc
cp calc_settings.json dist\PyCalc
cp -r SavedResults dist\PyCalc

