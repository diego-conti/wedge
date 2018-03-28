#ifndef DICTIONARY_PRODUCT_MANIFOLD_H
#define DICTIONARY_PRODUCT_MANIFOLD_H
using namespace Wedge;

class ProductManifold : public ConcreteManifold, public virtual Has_dTable {
	ExVector coordinates;
	ExVector CreateFrame(exvector oneforms, int no_coordinates)
	{
		ExVector frame;
  		frame.reserve(frame.size()+no_coordinates);
		frame=oneforms;
  		for (int i=1;i<=no_coordinates;i++)
			frame.push_back(DifferentialOneForm(Name("dx")(i)));
		return frame;
	}
	ExVector CreateFrame(exvector oneforms, ex coordinate)
	{
		ExVector frame;
		frame=oneforms;
		frame.push_back(DifferentialOneForm(Name("d"+ToString(coordinate))));
		LOG_INFO(frame);
		return frame;
	}
	void init (const Manifold& G)
	{
		VectorSpace<DifferentialForm> Lambda2G=G.pForms(2);
		Frame g(e().begin(),e().begin()+G.Dimension());
		VectorSpace<DifferentialForm> Lambda2GinProduct=::pForms(g,2);
		for (int i=1;i<=G.Dimension();++i) {
			ExVector comps=Lambda2G.Components(G.d(G.e(i)));
			ex dei;
			for (int n=1;n<=Lambda2G.Dimension();++n)
				dei+=comps(n)* Lambda2GinProduct.e(n);
			Has_dTable::Declare_d(e(i),dei);
		}
	}
public:
	ex x(OneBased i) {return coordinates(i);}

	ProductManifold(const Manifold& G, int extradimensions) : ConcreteManifold(CreateFrame(G.e(),extradimensions))
	{
		init(G);
		coordinates.reserve(extradimensions);
		for (int i=1;i<=extradimensions;i++)
		{
			Function x_i=Function(N.x(i));
			coordinates.push_back(x_i);
			Has_dTable::Declare_d(x_i,e(G.Dimension()+i));
			Has_dTable::Declare_d(e(G.Dimension()+i),0);
		}
	}
	ProductManifold(const Manifold& G, ex coordinate) : ConcreteManifold(CreateFrame(G.e(),coordinate))
	{
		init(G);
		coordinates.push_back(coordinate);
		Has_dTable::Declare_d(coordinate, e(G.Dimension()+1));
		Has_dTable::Declare_d(e(G.Dimension()+1),0);
	}

	ProductManifold(const Manifold& G, const Manifold& H) : ConcreteManifold(G.Dimension()+H.Dimension())
	{
		VectorSpace<DifferentialForm> Lambda2G=G.pForms(2), Lambda2H=H.pForms(2);
		Frame g(e().begin(),e().begin()+G.Dimension());
		Frame h(e().begin()+G.Dimension(),e().end());
		VectorSpace<DifferentialForm> Lambda2GinProduct=::pForms(g,2), Lambda2HinProduct=::pForms(h,2);
		for (int i=1;i<=G.Dimension();++i) {
			ExVector comps=Lambda2G.Components(G.d(G.e(i)));
			ex dei;
			for (int n=1;n<=Lambda2G.Dimension();++n)
				dei+=comps(n)* Lambda2GinProduct.e(n);
			Has_dTable::Declare_d(e(i),dei);
		}
		for (int i=1;i<=H.Dimension();++i) {
			ExVector comps=Lambda2H.Components(H.d(H.e(i)));
			ex dei;
			for (int n=1;n<=Lambda2H.Dimension();++n)
				dei+=comps(n)* Lambda2HinProduct.e(n);
			Has_dTable::Declare_d(e(i+G.Dimension()),dei);
		}
	}
};
#endif
