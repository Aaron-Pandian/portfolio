import csv
import PyPDF2
import re

# Task 1
data = open('find_the_link.csv',encoding='utf-8')
csv_data = csv.reader(data)
data_lines = list(csv_data)
google_link_unfinished = []
index_counter = 0

for row in data_lines:
    google_link_unfinished.append(row[index_counter])
    index_counter += 1
    
google_link = ''.join(google_link_unfinished)
# print(google_link)


# Task 2
f = open('Find_the_Phone_Number.pdf','rb')
pdf_reader = PyPDF2.PdfFileReader(f)

pages = []
for page in range(0,pdf_reader.getNumPages()):
    page = pdf_reader.getPage(page)
    page_text = page.extractText()
    pages.append(page_text)

phone_number = []
for page in pages:
    phone_find = re.search(r'\d{3}.\d{3}.\d{4}', page)

    try:
        phone_number.append(phone_find.group())
    except AttributeError:
        pass

# print(phone_number)

f.close()





