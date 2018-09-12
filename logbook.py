import pyqtgraph as pg
from subprocess import Popen, PIPE, check_output
from time import sleep
import os
from datetime import datetime
from xml.etree import ElementTree
from re import sub
from shutil import copy

# Shamelessly stolen from Shawn (thanks buddy). A lot of this is probably 
# unnecessary for my purposes but I'm too lazy to clean it up
def logbook(userText, titleText, textText, plotItem):
    curr_time = datetime.now()
    timeString = curr_time.strftime("%Y-%m-%dT%H:%M:%S")
    log_entry = ElementTree.Element(None)
    severity  = ElementTree.SubElement(log_entry, 'severity')
    location  = ElementTree.SubElement(log_entry, 'location')
    keywords  = ElementTree.SubElement(log_entry, 'keywords')
    time      = ElementTree.SubElement(log_entry, 'time')
    isodate   = ElementTree.SubElement(log_entry, 'isodate')
    log_user  = ElementTree.SubElement(log_entry, 'author')
    category  = ElementTree.SubElement(log_entry, 'category')
    title     = ElementTree.SubElement(log_entry, 'title')
    metainfo  = ElementTree.SubElement(log_entry, 'metainfo')
    imageFile = ElementTree.SubElement(log_entry, 'link')
    imageFile.text = timeString + '-00.ps'
    thumbnail = ElementTree.SubElement(log_entry, 'file')
    thumbnail.text = timeString + "-00.png"
    text      = ElementTree.SubElement(log_entry, 'text')
    log_entry.attrib['type'] = "LOGENTRY"
    category.text = "USERLOG"
    location.text = "not set"
    severity.text = "NONE"
    keywords.text = "none"
    time.text = curr_time.strftime("%H:%M:%S")
    isodate.text =  curr_time.strftime("%Y-%m-%d")
    metainfo.text = timeString + "-00.xml"
    fileName = "/tmp/" + metainfo.text
    fileName = fileName.rstrip(".xml")
    log_user.text = userText
    title.text = titleText
    text.text = textText
    
    # If field is truly empty, ElementTree leaves off tag entirely which 
    # causes logbook parser to fail
    if text.text == "": text.text = " "
     
    xmlFile = open(fileName+'.xml',"w")
    rawString = ElementTree.tostring(log_entry, 'utf-8')
    
    # Adds newline after each closing tag
    parsedString = sub(r'(?=<[^/].*>)','\n',rawString) 
    
    xmlString = parsedString[1:]
    xmlFile.write(xmlString)
    
    # Close with newline so cron job parses correctly
    xmlFile.write("\n")  
    
    xmlFile.close()
    exporter = pg.exporters.ImageExporter(plotItem)
    #exporter.parameters()['width'] = 550
    exporter.export(fileName+'.png')
    # PyQtGraph doesn't export PS files, so convert with linux 
    Popen('convert ' + fileName + '.png ' + fileName + '.ps', shell=True)
    
    sleep(0.3)
    Popen('convert ' + fileName + '.png -resize 500x500 ' + fileName
       + '.png', shell=True)
    sleep(0.35)
    path = "/u1/lcls/physics/logbook/data/"
    copy(fileName + '.ps', path)
    copy(fileName + '.png', path)
    copy(fileName + '.xml', path)

def MCCLog(tmpPNG, tmpPS, plotItem):
    exporter = pg.exporters.ImageExporter(plotItem)
    exporter.export(tmpPNG)
    Popen("convert "+ tmpPNG +" " + tmpPS, shell=True)
    sleep(0.1)
    printFile = "lpr -P" + 'elog_mcc' + " " + tmpPS
    os.system(printFile)
