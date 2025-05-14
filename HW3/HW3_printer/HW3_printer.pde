/* ||=========================================|| */
/* || PDA HW3 floorplan visualization program || */
/* ||=========================================|| */

void setup() {
    String caseName = "public1"; // change testcase
    float deadSpaceRatio = 0.15; // change dead space ratio
    
    int imageSize = 750;         // change image size
    size(750, 750);              // should be same as imageSize
    background(255);
    
    Printer printer = new Printer(imageSize, deadSpaceRatio);
    printer.readHardblocks("testcase/" + caseName + ".txt");
    printer.readFloorplan("output/" + caseName + ".out");
    printer.drawFloorplan("figure/" + caseName + ".png");
}

class Hardblock {
    String name;
    int w, h;
    int x, y;
    Boolean rotation;
    
    Hardblock(String name_, int w_, int h_) {
        name = name_;
        w = w_;
        h = h_;
        x = 0;
        y = 0;
        rotation = false;
    }
    
    void setPosition(int x_, int y_, Boolean rotation_) {
        x = x_;
        y = y_;
        rotation = rotation_;
    }
}

class Printer {
    int imageSize;
    float deadSpaceRatio;
    int totalBlockArea;
    ArrayList<Hardblock> hardblocks = new ArrayList<Hardblock>();
    
    IntDict strToIdx = new IntDict();
    
    Printer(int imageSize_, float deadSpaceRatio_) {
        imageSize = imageSize_;
        deadSpaceRatio = deadSpaceRatio_;
        totalBlockArea = 0;
    }
    
    void readHardblocks(String filename) {
        try {
            BufferedReader reader = createReader(filename);
            int idx = 0;
            String line;
            while ((line = reader.readLine()) != null) {
                String[] data = splitTokens(line, " ");
                if (data.length > 0 && data[0].equals("HardBlock")) {
                    String name = data[1];
                    int w = int(data[2]);
                    int h = int(data[3]);
                    totalBlockArea += w * h;
                    hardblocks.add(new Hardblock(name, w, h));
                    strToIdx.set(name, idx++);
                }
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    void readFloorplan(String filename) {
        try {
            BufferedReader reader = createReader(filename);
            String line;
            while ((line = reader.readLine()) != null) {
                String[] data = splitTokens(line, " ");
                if (data.length == 4) {
                    String name = data[0];
                    int x = int(data[1]);
                    int y = int(data[2]);
                    boolean rotation = int(data[3]) == 1;
                    if (strToIdx.hasKey(name)) {
                        hardblocks.get(strToIdx.get(name)).setPosition(x, y, rotation);
                    }
                }
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    
    int getOutlineSize() {
        return floor(sqrt(totalBlockArea * (1 + deadSpaceRatio)));
    }
    
    int calX(int x) {
        return x + 1;
    }
    
    int calY(int y, int h) {
        return imageSize - y - h - 2;
    }
    
    void drawOutline() {
        int size = getOutlineSize();
        stroke(255, 0, 0);
        noFill();
        rect(calX(0), calY(0, size), size, size);
    }
    
    void drawHardblocks() {
        for (Hardblock it : hardblocks) {
            int w, h, x, y;
            if (it.rotation) {
                w = it.h;
                h = it.w;
                x = calX(it.x);
                y = calY(it.y, it.w);
            }
            else {
                w = it.w;
                h = it.h;
                x = calX(it.x);
                y = calY(it.y, it.h);
            }
            
            stroke(0);
            fill(200);
            rect(x, y, w, h);
            
            fill(0);
            textSize(8);
            text(it.name, x + 2, y + 9);
        }
    }
    
    void drawFloorplan(String filename) {
        drawHardblocks();
        drawOutline();
        save(filename);
    }
}
