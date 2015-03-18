/*
 * Window.cpp
 *
 *  Created on: Nov 26, 2012
 *      Author: jagiella
 */




#include "Window.hpp"

#include <QtGui/QPainter>
#include <QtGui/QMouseEvent>
#include <QtCore/QTimer>

#include <math.h>

Window::Window( Cells *cells, Boxes<Agent*> *box, float *time,  QWidget *parent ) : QWidget( parent )
{
  r = 0;

  zoom=1;
  this->setWindowTitle( "[Zoom: use mouse wheel]");

  timer = new QTimer( this );
  timer->setInterval( 50 );

  connect( timer, SIGNAL(timeout()), this, SLOT(timeout()) );

  this->cells = cells;
  this->box = box;
  this->time = time;

  timer->start();
}

Window::Window( AgentList<> *cells, Boxes<Agent*> *box, float *time,  QWidget *parent ) : QWidget( parent )
{
	printf("Create Window\n");
  r = 0;

  zoom=0.75;
  this->setWindowTitle( "[Zoom: use mouse wheel]");

  timer = new QTimer( this );
  timer->setInterval( 50 );

  connect( timer, SIGNAL(timeout()), this, SLOT(timeout()) );

  this->cells = 0;
  this->al = cells;
  this->box = box;
  this->time = time;

  image_idx = -1;

  timer->start();
}

QSize Window::sizeHint() const
{
  return QSize( 400, 400 );
}

void Window::timeout()
{
  if( r == 0 )
  {
    x = mx;
    y = my;

    color = QColor( qrand()%256, qrand()%256, qrand()%256 );
    growing = QColor( 255, 0, 0 );
    quiescent = QColor( 0, 255,0 );
  }

  int dx = mx-x;
  int dy = my-y;

  if( dx*dx+dy*dy <= r*r )
    r++;
  else
    r--;

  if(image_idx < (int)(*time*2.)){
	  char filename[512];
	  sprintf(filename, "test%i.png", (int)(*time*2.));
	  saveToFile(filename);
	  image_idx++;
  }
  update();
}

#define MIN(a,b) (a<b?a:b)

void Window::paintEvent( QPaintEvent* )
{
	//printf("Paint Window\n");
  if( r > 0 )
  {
    QPainter painter( this );

    painter.setRenderHint( QPainter::Antialiasing );

    QColor border( 0,0,0);
    painter.setPen( border );


    int N;
    if(cells)
    	N=cells->N;
    else
    	N=al->size();
    float scale = zoom * MIN(this->width(),this->height()) / 15. / sqrt(N*4/3.142);

    if(cells)
    for( int i=0; i<cells->N; i++){
    	//QColor *cellcolor;
    	if( cells->type[i] == Cells::barrier){
    		color.setRgb( 0,0,255);
    		painter.setBrush( color );
    	}
    	else
    	if( cells->cycleTime[i]>0 && !cells->overlapp[i]){
    		//painter.setBrush( growing );
    		color.setRgb(MIN( cells->cycleTime[i]/24.*255.,255),0,0);
    		painter.setBrush( color );
    		//cellcolor = new QColor( 1,0,0);
    	}
    	else
    		painter.setBrush( quiescent );
    		//color.setRgb(0,1,0);
    		//cellcolor = new QColor( 0,1,0);
    	//painter.setBrush( color );
    	painter.drawEllipse(
    			scale*(cells->position[i][0]-cells->radius[i]-box->X*box->latticeConstant/2)  + this->width()/2,
    			scale*(cells->position[i][1]-cells->radius[i]-box->Y*box->latticeConstant/2)  + this->height()/2,
    			scale*(cells->radius[i]*2), scale*(cells->radius[i]*2) );
    	//delete( cellcolor);
    }
    else
        for( int i=0; i<al->size(); i++){
        	if( al->at(i)->type() == 2){
        		color.setRgb( 0,0,255);
        		painter.setBrush( color );
        	}else if( al->at(i)->pressure()<4e2)
        		painter.setBrush( growing );
        	else
        		painter.setBrush( quiescent );
        	painter.drawEllipse(
        			scale*(al->at(i)->position()[0] - al->at(i)->radius() - box->X*box->latticeConstant/2)  + this->width()/2,
        			scale*(al->at(i)->position()[1] - al->at(i)->radius() - box->Y*box->latticeConstant/2)  + this->height()/2,
        			scale*(al->at(i)->radius()*2), scale*(al->at(i)->radius()*2) );
        }

    char titel[512];
	sprintf(titel, "%.0lf days, %.0lf hours", floor(*time/24.), *time-floor(*time/24.)*24.);
    painter.drawText( QPoint(10,20), titel );

  }

}

void Window::mousePressEvent( QMouseEvent *e )
{
  mx = e->x();
  my = e->y();

  if( timer->isActive()){
	  timer->stop();
		char titel[512];
		sprintf(titel, "[Click to Resume]");
		this->setWindowTitle( titel);

  }
  else{
	  timer->start();
		char titel[512];
		sprintf(titel, "[Zoom: %.0f %%]", zoom*100);
		this->setWindowTitle( titel);

  }
}

void Window::mouseMoveEvent( QMouseEvent *e )
{
  mx = e->x();
  my = e->y();
}

void Window::mouseReleaseEvent( QMouseEvent *e )
{
  //timer->stop();
	UNUSED(e);
}

void Window::wheelEvent ( QWheelEvent * event )
{
	//fprintf(stderr, "\nzoom ->  %i (%f)\n", zoom, event->delta()/500.);

	zoom += zoom * event->delta()/500.;

	char titel[512];
	sprintf(titel, "[Zoom: %.0f %%]", zoom*100);
	this->setWindowTitle( titel);
}

void Window::saveToFile ( QString saveFilename ){
	QString saveExtension = "PNG";
	int pos = saveFilename.lastIndexOf('.');
	if (pos >= 0)
	    saveExtension = saveFilename.mid(pos + 1);

	if(!QPixmap::grabWidget( this).save(saveFilename, qPrintable(saveExtension)))
	{
	    // since you have a widget, just use grabWidget() here. winId() would possibly have
	    // portability issues on other platforms.  qPrintable(saveExtension) is effectively
	    // the same as saveExtension.toLocal8Bit().constData()

	    //QMessageBox::warning(this, "File could not be saved", "ok", QMessageBox::Ok);
	}
}
