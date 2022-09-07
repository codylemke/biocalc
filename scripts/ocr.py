"""Returns text from an input PDF

Descripion

Usage Example:
    example
"""
# Standard Libraries
import pathlib
# Third Party Libraries
import pdf2image
import PIL
import pytesseract
# Global Constants
ROOT_DIR = pathlib.Path(__file__).resolve().parent


# SCRIPT ----------------------------------------------------------------------
if __name__ == '__main__':
    pdf_file = ROOT_DIR / 'inputs' / 'input.pdf'
    print('starting conversion')
    images = pdf2image.convert_from_path(pdf_file)
    print('starting loop')
    text = list()
    for iteration, image in enumerate(images):
        print(f'iteration: {iteration}')
        text.append(pytesseract.image_to_string(image))
        # print(text)
        # image.save(f'page{iteration}.jpg', 'JPEG')
    # text = pytesseract.image_to_string(PIL.Image.open(pdf_file))
    # print(text)
    text_file = ROOT_DIR / 'outputs' / 'pdf.txt'
    text_file.touch()
    with open(text_file, 'w') as file:
        file.write(''.join(text))