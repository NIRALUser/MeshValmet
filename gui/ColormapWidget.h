#ifndef COLORMAPWIDGET_H
#define COLORMAPWIDGET_H

#include <qwidget.h>

class ColormapWidget : public QWidget
{
Q_OBJECT
public:
  ColormapWidget( QWidget *parent=0, const char *name=0 );
  ~ColormapWidget();
  
protected:
  void paintEvent(QPaintEvent *); 
  
private:

struct look
    {
      float R;
      float G;
      float B;
    };
look looktable[256];

};

#endif
