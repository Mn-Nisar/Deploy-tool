import io
import re
import PIL
from PIL import ImageDraw, ImageColor
import base64
import requests
from io import BytesIO

def parse_conf(conf_result):
    for line in conf_result.splitlines():
        linelist = line.strip().split('\t')
        res = [each for each in re.split(
            r'[\s,\(\)]', linelist[0]) if each]
        shape = res[0]
        position = list(map(int, res[1:]))

        url = linelist[1]
        title = linelist[2]
        yield shape, position, url, title


def get_gene_color(title, gene_color, conflict_color='green'):
    color = None
    result = re.findall(r'([^\s]+) \((.+?)\)', title)
    if result:
        gene_in_title = [each for part in result for each in part]
        for gene in gene_in_title:
            if gene_color.get(gene):
                # conflict color in a gene family
                if color and color != gene_color[gene]:
                    color = conflict_color
                else:
                    color = gene_color[gene]
    return color


def build_png(img , conf_data, gene_color):
        im = PIL.Image.open(img)
        draw = ImageDraw.Draw(im)
        for shape, position, _, title in conf_data:
            color = get_gene_color(title, gene_color)
            if not color:
                continue

            try:
                color_rgba = ImageColor.getcolor(color, 'RGBA')
            except:
                color_rgba = ImageColor.getcolor(f'#{color}', 'RGBA')

            if not color_rgba:
                continue

            if shape == 'rect':
                X, Y, RX, RY = position
                for x in range(X, RX):
                    for y in range(Y, RY):
                        # pixel > 0 means this point is not black
                        if im.getpixel((x, y))[0] > 0:
                            ImageDraw.floodfill(
                                im, xy=(x, y), value=color_rgba)
            elif shape =="circ":
                CX, CY, R = position
                X, Y, RX, RY = CX - R, CY - R, CX + R, CY + R
                draw.ellipse([(X, Y), (RX, RY)], fill=color_rgba)

        file = io.BytesIO()
        im.save(file, format='png')
        # im.save('mypathh.png')
        file.seek(0)
        image_png = file.getvalue()
        pathway = base64.b64encode(image_png)
        pathway = pathway.decode('utf-8')
        file.close()
        return pathway


def draw_pathway(req_pathway , color_dict , request):
    req_pathway = req_pathway.split(':')[-1]
    req_pathway = req_pathway.strip()
    req_pathway = 'map'+req_pathway
    image_link = "https://rest.kegg.jp/get/"+req_pathway+"/image"
    conf_link = "https://rest.kegg.jp/get/"+req_pathway+"/conf"

    image_file = requests.get(image_link, params=request.GET)

    conf_file = requests.get(conf_link, params=request.GET)

    if image_file.status_code == 200 and conf_file.status_code == 200:

        img = BytesIO(image_file.content)

        conf_data = list(parse_conf(conf_file.text))
        final_image = build_png(img , conf_data,  color_dict)

        return final_image
