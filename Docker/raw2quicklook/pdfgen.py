import os
import math
import io
from io import BytesIO
from PIL import Image
from pathlib import Path

import PyPDF2
from reportlab.lib import colors
from reportlab.pdfbase.pdfmetrics import stringWidth
from reportlab.pdfgen.canvas import Canvas
from reportlab.platypus import Table, TableStyle

class pdfGen: 

    def __init__(self,title,folder,imagesPerPage=3,columns=1):
        self.title=title
        self.path = os.path.join(folder, Path(title).stem)
        self.path = folder
        self.constructor()
        self.pdfFile = PyPDF2.PdfWriter()
        self.columns=columns
        self.imagesPerPage=imagesPerPage
        self.rows = math.ceil(self.imagesPerPage / self.columns)
        self.imgCount=0
        self.pages=[]
        self.pdfWidth = 2000
        self.pdfHeight = 1000

    def constructor(self): #TODO 
        os.makedirs(self.path,exist_ok=True)
        # isExist = os.path.exists(self.path)
        # if not isExist:
        #     os.mkdir(self.path)
        #     os.mkdir(self.path + '/plots')
        # else: 
        #     os.makedirs(self.path + '/plots',exist_ok=True)

    def addPlot(self,fig):
        # fig.savefig(self.path+'/plots/'+name)
        # fig.write_image(self.path+'/plots/'+name, engine='kaleido', scale=2, format='png', width=1200, height=600)
        # image = Image.open(self.path+'/plots/'+name)
        # image = image.convert('RGB')
        imgData = BytesIO()
        fig.update_layout(height=600, width=1200)
        fig.write_image(imgData, engine='kaleido', scale=2, format='png')
        image = Image.open(imgData)
        image = image.convert('RGB')
        pdf_output = BytesIO()
        image.save(pdf_output, format='PDF')
        pdf = PyPDF2.PdfReader(pdf_output)
        self.addImg(pdf.pages[0])

    def addImg(self,image,xCoor=float('NaN'),yCoor=float('NaN'),pageNumber=float('NaN')):
        imageBox = image.mediabox
        if(self.imgCount%self.imagesPerPage==0 and len(self.pages)<=self.imgCount//self.imagesPerPage):
            self.pdfWidth = self.columns * imageBox.width
            self.pdfHeight = self.rows * imageBox.height
            self.addBlankPage()
        if(math.isnan(xCoor) or math.isnan(yCoor) or math.isnan(pageNumber)):
            pageBox = self.pages[self.imgCount//self.imagesPerPage].mediabox
            tempPage = PyPDF2.PdfWriter().add_blank_page(pageBox.width, pageBox.height)
            tempPage.merge_page(image)
            tempPage.add_transformation(PyPDF2.Transformation().translate(((self.imgCount%self.imagesPerPage) % self.columns) * imageBox.width, (self.rows - 1 - (self.imgCount%self.imagesPerPage) // self.columns) * imageBox.height))
            self.pages[self.imgCount//self.imagesPerPage].merge_page(tempPage)
            self.imgCount+=1
        elif(len(self.pages)>pageNumber):
            pageBox = self.pages[int(pageNumber)].mediabox
            tempPage = PyPDF2.PdfWriter().add_blank_page(pageBox.width, pageBox.height)
            tempPage.merge_page(image)
            tempPage.add_transformation(PyPDF2.Transformation().translate(xCoor,pageBox.height-yCoor))
            self.pages[int(pageNumber)].merge_page(tempPage)
        else:
            print("Cannot add the image, wrong page number")

    def addBlankPage(self):
        self.pages.append(PyPDF2.PdfWriter().add_blank_page(self.pdfWidth, self.pdfHeight))

    def addText(self,text,xCoor,yCoor,pageNumber=float('NaN'),fontSize=60,alignment='none'):
        if(math.isnan(pageNumber)):
            pageNumber=len(self.pages)-1
        packet = io.BytesIO()
        c = Canvas(packet, pagesize=(self.pages[pageNumber].mediabox.width, self.pages[pageNumber].mediabox.height))
        textobject = c.beginText()
        textobject.setFont('Helvetica', fontSize)
        if(alignment=='center'):
            textobject.setTextOrigin((int(self.pages[pageNumber].mediabox.width/2)-int(stringWidth(text,'Helvetica', fontSize)/2)), int(self.pdfHeight-yCoor))
        elif(alignment=='right'):
            textobject.setTextOrigin((int(self.pages[pageNumber].mediabox.width)-int(stringWidth(text,'Helvetica', fontSize))), int(self.pdfHeight-yCoor))
        elif(alignment=='left'):
            textobject.setTextOrigin(0, int(self.pdfHeight - yCoor))
        else:
            textobject.setTextOrigin(xCoor, int(self.pdfHeight-yCoor))
        textobject.textLines(text)
        c.drawText(textobject)
        c.save()
        packet.seek(0)
        pagePdf = PyPDF2.PdfReader(packet)
        page = pagePdf.pages[0]
        self.pages[pageNumber].merge_page(page)

    def addTable(self,data,xCoor,yCoor,pageNumber=float('NaN')):
        if(math.isnan(pageNumber)):
            pageNumber=len(self.pages)-1
        packet = io.BytesIO()
        c = Canvas(packet, pagesize=(self.pages[pageNumber].mediabox.width, self.pages[pageNumber].mediabox.height))
        table=Table(data)
        table.setStyle(TableStyle([
                                ('BACKGROUND', (0, 0), (-1, 0),colors.powderblue),
                                ('BACKGROUND', (0, 1), (0, -1),colors.cadetblue),
                                ('LEADING', (0, 0), (-1, -1), 80),
                                ('FONTSIZE', (0, 0), (-1, -1), 40),
                                ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
                                ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
                                ('ALIGN', (0,0), (-1,-1), 'CENTER'),
                                ('BOX', (0, 0), (-1, -1), 0.25, colors.black),
        ]))
        table.wrapOn(c, 500, 500)
        table.drawOn(c, xCoor, int(self.pdfHeight-yCoor))
        c.save()
        packet.seek(0)
        pagePdf = PyPDF2.PdfReader(packet)
        page = pagePdf.pages[0]
        self.pages[pageNumber].merge_page(page)

    def pdfFileGenerator(self):
        for i,item in enumerate(self.pages):
            self.addText(Path(self.title).stem,0,int(item.mediabox.height)-10,i,fontSize=20,alignment='right')
            self.pdfFile.add_page(item)

        with open(self.path + '/' + Path(self.title).stem + '.pdf', 'wb') as f:
            self.pdfFile.write(f)
    