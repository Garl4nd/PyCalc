echo y | pyinstaller --onedir --icon="pi.ico" PyCalc.pyw
xcopy pi.ico dist\PyCalc
